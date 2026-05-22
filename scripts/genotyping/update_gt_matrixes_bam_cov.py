import pandas as pd
import sys
import os
import numpy as np


GROUP2_SAMPLES = {'RCC1749', 'RCC3052'}


gt_matrix_file = sys.argv[1]
bed_path = sys.argv[2]
gt_matrix_out = sys.argv[3]

gt_df = pd.read_csv(gt_matrix_file, sep="\t", header=0)


def _clean_status(value):
    if pd.isna(value) or value == '' or value == 'NA':
        return 'NA'
    return str(value)


def _group_pattern(present, absent, missing):
    callable_count = present + absent
    if callable_count == 0:
        return 'no_call'
    if present == 0:
        return 'absent'
    if absent == 0 and missing == 0:
        return 'fixed_present'
    if absent == 0:
        return 'present_with_missing'
    if present == 1:
        return 'singleton_present'
    if absent == 1:
        return 'singleton_absent'
    return 'polymorphic'


def _within_orthology_confidence(status):
    status = _clean_status(status)
    if status in {'consistent', 'singleton'}:
        return 'high'
    if status == 'low_identity':
        return 'low_identity'
    if status in {'discordant', 'uncertain'}:
        return status
    return 'unknown'


def _cross_group_origin(status):
    status = _clean_status(status)
    if status in {'ancestral', 'likely_ancestral', 'ancestral_low_identity'}:
        return 'ancestral'
    if status in {'independent', 'likely_independent'}:
        return 'independent'
    if status == 'uncertain':
        return 'ambiguous'
    return 'not_applicable'


def _cross_group_confidence(status):
    status = _clean_status(status)
    if status in {'ancestral', 'independent'}:
        return 'high'
    if status in {'likely_ancestral', 'likely_independent'}:
        return 'moderate'
    if status == 'ancestral_low_identity':
        return 'low_identity'
    if status == 'uncertain':
        return 'low'
    return 'not_applicable'


def add_derived_classification_columns(df):
    """Add explicit, paper-facing classification columns.

    The legacy status columns are intentionally preserved for downstream
    compatibility. These derived columns separate per-group callability and
    frequency pattern from orthology confidence and cross-group origin.
    """
    group2 = set(GROUP2_SAMPLES)
    all_samples = sorted(df['sample'].dropna().unique())
    group1 = [s for s in all_samples if s not in group2]
    group2 = [s for s in all_samples if s in group2]

    derived = {}
    first_rows = df.drop_duplicates('ortholog_id').set_index('ortholog_id')

    for oid, group in df.groupby('ortholog_id', sort=False):
        record = {}
        for label, samples in [('group1', group1), ('group2', group2)]:
            calls = group.loc[group['sample'].isin(samples), 'presence']
            present = int((calls == 1).sum())
            absent = int((calls == 2).sum())
            missing = int((calls == 3).sum())
            callable_count = present + absent
            n_samples = len(calls)
            record[f'{label}_present_count'] = present
            record[f'{label}_absent_count'] = absent
            record[f'{label}_missing_count'] = missing
            record[f'{label}_callable_count'] = callable_count
            record[f'{label}_n_samples'] = n_samples
            record[f'{label}_callability'] = (
                'complete' if missing == 0 else
                'partial' if callable_count > 0 else
                'none'
            )
            record[f'{label}_pattern'] = _group_pattern(present, absent, missing)

        first = first_rows.loc[oid]
        within_status = first.get('within_group_status', '')
        cross_status = first.get('cross_group_status', '')
        record['within_group_orthology_confidence'] = (
            _within_orthology_confidence(within_status)
        )
        record['cross_group_origin'] = _cross_group_origin(cross_status)
        record['cross_group_confidence'] = _cross_group_confidence(cross_status)
        record['present_in_both_groups'] = (
            record['group1_present_count'] > 0 and
            record['group2_present_count'] > 0
        )
        derived[oid] = record

    derived_df = pd.DataFrame.from_dict(derived, orient='index')
    derived_df.index.name = 'ortholog_id'

    # Replace these derived columns if re-running on an already annotated file.
    for col in derived_df.columns:
        if col in df.columns:
            df = df.drop(columns=[col])

    return df.merge(derived_df.reset_index(), on='ortholog_id', how='left')

# update genotype matrix with new calls
for file in os.listdir(bed_path):
    if file.endswith("coverage_calls.bed"):
        print("Reading...", file)
        file_path = os.path.join(bed_path, file)
        sample = file.split(".")[0]
        bed_df = pd.read_csv(file_path, sep="\t", header=0)

        for index, row in bed_df.iterrows():
            chrom, start, end, ortho_id, new_call = row['contig'], row['start'], row['end'], row['ortholog_id'], int(row['new_call'])
            # update call df with new call
            gt_df.loc[(gt_df['contig'] == chrom) & (gt_df['start'] == start) & (gt_df['end'] == end) & (gt_df['ortholog_id'] == ortho_id), 'presence'] = new_call


# Reconcile classification labels with post-coverage presence state.
# Coverage can flip 1 -> 2/3, which can invalidate the pre-coverage labels:
#   - A multi-member "consistent" group with only 1 remaining member is now
#     a "singleton"; identity columns lose meaning
#   - A cross-group (G1+G2) ortholog that lost all members of one clade is
#     now within-group-only, so cross_group_status should be NA
#   - A group with zero presence=1 members is uninformative and gets dropped
# Multi-member groups that lost some but kept >=2 members stay as-is — their
# labels remain biologically accurate for the remaining present set.
if 'within_group_status' in gt_df.columns:
    g2_mask = gt_df['sample'].isin(GROUP2_SAMPLES)
    present_mask = gt_df['presence'] == 1

    present_rows = gt_df[present_mask]
    n_present = present_rows.groupby('ortholog_id').size()
    n_g1_present = present_rows[~present_rows['sample'].isin(GROUP2_SAMPLES)] \
        .groupby('ortholog_id').size()
    n_g2_present = present_rows[present_rows['sample'].isin(GROUP2_SAMPLES)] \
        .groupby('ortholog_id').size()

    all_oids = set(gt_df['ortholog_id'].unique())
    empty_groups = all_oids - set(n_present.index)
    singleton_oids = set(n_present[n_present == 1].index)

    # Groups whose pre-coverage cross_group_status was non-NA but which now
    # lack members in one clade: downgrade cross_group_status to NA.
    def _is_na_cross(v):
        return pd.isna(v) or v == '' or v == 'NA'

    first_rows = gt_df.drop_duplicates('ortholog_id').set_index('ortholog_id')
    cross_to_na_oids = set()
    if 'cross_group_status' in gt_df.columns:
        for oid, row in first_rows.iterrows():
            if oid in empty_groups:
                continue
            if _is_na_cross(row.get('cross_group_status', '')):
                continue
            has_g1 = n_g1_present.get(oid, 0) > 0
            has_g2 = n_g2_present.get(oid, 0) > 0
            if not (has_g1 and has_g2):
                cross_to_na_oids.add(oid)

    # Identity columns load as float64; cast to object so we can write
    # empty strings for cleared values
    for col in ('within_group_identity', 'cross_group_identity'):
        if col in gt_df.columns:
            gt_df[col] = gt_df[col].astype(object)

    # Apply: within-group-status = singleton (and clear within_group_identity)
    if singleton_oids:
        mask = (gt_df['ortholog_id'].isin(singleton_oids) &
                (gt_df['within_group_status'] != 'singleton'))
        n_downgraded = gt_df.loc[mask, 'ortholog_id'].nunique()
        gt_df.loc[mask, 'within_group_status'] = 'singleton'
        if 'within_group_identity' in gt_df.columns:
            gt_df.loc[gt_df['ortholog_id'].isin(singleton_oids),
                      'within_group_identity'] = ''
        print(f"Downgraded {n_downgraded} groups to singleton "
              f"(after coverage left them with 1 present member)")

    # Apply: cross_group_status = NA (and clear cross_group_identity)
    if cross_to_na_oids:
        mask = gt_df['ortholog_id'].isin(cross_to_na_oids)
        gt_df.loc[mask, 'cross_group_status'] = ''
        if 'cross_group_identity' in gt_df.columns:
            gt_df.loc[mask, 'cross_group_identity'] = ''
        print(f"Reset cross_group_status to NA for {len(cross_to_na_oids)} "
              f"groups (one clade fully absent after coverage)")

    # Drop groups with zero presence=1 members after coverage
    if empty_groups:
        n_empty_rows = gt_df['ortholog_id'].isin(empty_groups).sum()
        gt_df = gt_df[~gt_df['ortholog_id'].isin(empty_groups)].copy()
        print(f"Dropped {len(empty_groups)} empty groups "
              f"({n_empty_rows} rows) with zero present members after coverage")

# Convert numeric columns to integer type, handling NaN values
if 'start' in gt_df.columns:
    # Replace NaN values with a placeholder (-1)
    gt_df['start'] = gt_df['start'].fillna(-1).astype(int)
    # If you prefer NaN values to remain as NaN, use Int64 type instead:
    # gt_df['start'] = gt_df['start'].astype('Int64')

if 'end' in gt_df.columns:
    gt_df['end'] = gt_df['end'].fillna(-1).astype(int)
    # Alternative: gt_df['end'] = gt_df['end'].astype('Int64')

if 'family' in gt_df.columns:
    # First check if column is numeric
    if pd.api.types.is_numeric_dtype(gt_df['family']):
        gt_df['family'] = gt_df['family'].fillna(-1).astype(int)
        # Alternative: gt_df['family'] = gt_df['family'].astype('Int64')
    else:
        # For non-numeric family column, convert only where possible
        try:
            # For string columns with numeric content
            numeric_mask = pd.to_numeric(gt_df['family'], errors='coerce').notna()
            gt_df.loc[numeric_mask, 'family'] = pd.to_numeric(gt_df.loc[numeric_mask, 'family']).astype(int)
        except Exception as e:
            print(f"Warning: Could not fully convert family column: {e}")
            # Keep as is if conversion fails

gt_df = add_derived_classification_columns(gt_df)

# Debug information
print(f"Column dtypes after conversion: {gt_df.dtypes}")
print(f"NaN values in start: {gt_df['start'].isna().sum()}")
print(f"NaN values in end: {gt_df['end'].isna().sum()}")
if 'family' in gt_df.columns:
    print(f"NaN values in family: {gt_df['family'].isna().sum()}")

# Save the updated dataframe
gt_df.to_csv(gt_matrix_out, sep="\t", header=True, index=False)
print(f"Updated matrix saved to {gt_matrix_out}")
