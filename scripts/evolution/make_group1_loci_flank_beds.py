#!/usr/bin/env python3

import sys
import os
import pandas as pd

def main():
    if len(sys.argv) != 6:
        print("Usage: python make_group1_loci_flank_beds.py <metadata_file> <sample_name> <flank_length> <left_bed> <right_bed>")
        sys.exit(1)
    
    metadata_file = sys.argv[1]
    sample_name = sys.argv[2]
    flank_length = int(sys.argv[3])
    left_bed = sys.argv[4]
    right_bed = sys.argv[5]

    # Read genotype matrix
    meta_df = pd.read_csv(metadata_file, sep="\t", header=0)

    # Filter for the specific sample
    sample_df = meta_df[meta_df['sample'] == sample_name]
    sample_df = sample_df[sample_df['presence'] != 3]  # Remove missing data
    sample_df = sample_df[sample_df['start'] >= flank_length]  # Ensure sufficient flanking space
    sample_df['start'] = sample_df['start'].astype(int)
    sample_df['end'] = sample_df['end'].astype(int)
    sample_df = sample_df.sort_values(by=['contig', 'start'])

    with open(left_bed, 'w') as l, open(right_bed, 'w') as r:
        for index, row in sample_df.iterrows():
            contig = row['contig']
            
            # Handle sequences that are reverse complemented
            # Maintain biological directionality (5' and 3' flanks)
            if row['orientation'] == 'reverse':
                # For reverse orientation: 
                # - biological left (5') flank extends upstream from sequence end
                # - biological right (3') flank extends downstream from sequence start
                left_start = row['end']
                left_end = row['end'] + flank_length
                right_start = row['start'] - flank_length
                right_end = row['start']
            else:
                # For forward orientation:
                # - biological left (5') flank extends upstream from sequence start  
                # - biological right (3') flank extends downstream from sequence end
                left_start = row['start'] - flank_length
                left_end = row['start']
                right_start = row['end']
                right_end = row['end'] + flank_length
                
            ortholog_id = row['ortholog_id']

            l.write(f"{contig}\t{left_start}\t{left_end}\t{ortholog_id}\n")
            r.write(f"{contig}\t{right_start}\t{right_end}\t{ortholog_id}\n")

if __name__ == "__main__":
    main()