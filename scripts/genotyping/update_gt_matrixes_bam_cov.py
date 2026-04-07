import pandas as pd
import sys
import os
import numpy as np


gt_matrix_file = sys.argv[1]
bed_path = sys.argv[2]
gt_matrix_out = sys.argv[3]

gt_df = pd.read_csv(gt_matrix_file, sep="\t", header=0)

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

# Debug information
print(f"Column dtypes after conversion: {gt_df.dtypes}")
print(f"NaN values in start: {gt_df['start'].isna().sum()}")
print(f"NaN values in end: {gt_df['end'].isna().sum()}")
if 'family' in gt_df.columns:
    print(f"NaN values in family: {gt_df['family'].isna().sum()}")

# Save the updated dataframe
gt_df.to_csv(gt_matrix_out, sep="\t", header=True, index=False)
print(f"Updated matrix saved to {gt_matrix_out}")
