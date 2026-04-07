#!/usr/bin/env python3
import pandas as pd
import sys

bed_file = sys.argv[1]
filtered_bed_file = sys.argv[2]

# Define column names 
columns = [
    'sseqid', 'start', 'end', 'qseqid', 'pident', 'qlen', 'slen',
    'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'evalue'
]

# Read the input file
df = pd.read_csv(bed_file, sep='\t', names=columns)

# -------------------------------------------------------------------------
# 1) Vectorized filtering to handle the "pident >= 90" logic

# Calculate the max pident for each (sseqid, start, end)
df['pident_max'] = df.groupby(['sseqid', 'start', 'end'])['pident'].transform('max')

# Create a mask:
#   - if max pident >= 90, keep only rows with pident >= 90
#   - otherwise, keep only rows with pident == pident_max
mask = (
    ((df['pident_max'] >= 90) & (df['pident'] >= 90)) |
    ((df['pident_max'] < 90) & (df['pident'] == df['pident_max']))
)
filtered_df = df[mask].copy()

# We no longer need the helper column
filtered_df.drop(columns='pident_max', inplace=True)

# -------------------------------------------------------------------------
# 2) Single-pass overlap filter

# First, sort by sseqid, start, end
filtered_df.sort_values(by=['sseqid', 'start', 'end'], inplace=True)

to_keep = []
current_sseqid = None
best_row = None
previous_end = -1

# Iterate over rows in sorted order
for idx, row in filtered_df.iterrows():
    # If we are switching to a new sseqid, finalize the previous best
    if row['sseqid'] != current_sseqid:
        if best_row is not None:
            to_keep.append(best_row.name)
        current_sseqid = row['sseqid']
        best_row = row
        previous_end = row['end']
        continue

    # Same sseqid: check overlap
    if row['start'] <= previous_end:
        # Overlap. Keep the row with the longest qlen, then highest pident
        if (best_row is None 
            or (row['qlen'] > best_row['qlen']) 
            or (row['qlen'] == best_row['qlen'] and row['pident'] > best_row['pident'])):
            best_row = row
    else:
        # No overlap: finalize the previous best
        to_keep.append(best_row.name)
        best_row = row

    # Update previous_end
    previous_end = max(previous_end, row['end'])

# Finalize the very last best row
if best_row is not None:
    to_keep.append(best_row.name)

# Extract the rows we decided to keep
final_df = filtered_df.loc[to_keep]

# Write the final results
final_df.to_csv(filtered_bed_file, sep='\t', index=False, header=False)