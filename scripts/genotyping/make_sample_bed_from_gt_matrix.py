import pandas as pd 
import sys

gt_matrix = sys.argv[1]
sample = sys.argv[2]
bed = sys.argv[3]

df = pd.read_csv(gt_matrix, sep="\t", header=0)
df = df[df['sample'] == sample]
df = df[df['presence'] == 1]

#reorder columns
df = df[['contig', 'start', 'end', 'sample', 'ortholog_id', 'sequence_id', 'family', 'gene', 'splice_site', 'presence', 'orientation', 'left_reverse', 'right_reverse']]
df = df.sort_values(by=['contig', 'start'])
cols  = ['start', 'end', 'family']
df[cols] = df[cols].astype(int)

df.to_csv(bed, sep="\t", header=True, index=False)