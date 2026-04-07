from Bio import SeqIO
import sys

def extract_flanks(input_fa, left_out, right_out):
    with open(left_out, 'w') as left_f, open(right_out, 'w') as right_f:
        for record in SeqIO.parse(input_fa, 'fasta'):
            seq = str(record.seq)
            left_flank = seq[:100]
            right_flank = seq[-100:]
            
            left_f.write(f">{record.id}_left\n{left_flank}\n")
            right_f.write(f">{record.id}_right\n{right_flank}\n")

if __name__ == "__main__":
    extract_flanks(sys.argv[1], sys.argv[2], sys.argv[3])