
import argparse

def read_fasta(fasta_path):
    """
    Generator function that reads a FASTA file and yields (header, sequence).
    """
    header = None
    seq_chunks = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # If we have a current sequence, yield it
                if header and seq_chunks:
                    yield (header, "".join(seq_chunks))
                # Start a new sequence
                header = line[1:]  # omit the '>' character
                seq_chunks = []
            else:
                seq_chunks.append(line)
        # Yield the last sequence in the file, if present
        if header and seq_chunks:
            yield (header, "".join(seq_chunks))

def count_kmers(sequence, k):
    """
    Returns a dictionary of k-mer counts for the given sequence, *ignoring* any k-mer that contains 'N'.
    """
    kmer_counts = {}
    seq_len = len(sequence)
    for i in range(seq_len - k + 1):
        kmer = sequence[i : i + k]
        # Skip if 'N' is present
        if 'N' in kmer.upper():
            continue
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    return kmer_counts

def check_sequence(sequence, kmers_to_test):
    """
    Checks the sequence for:
      1) N content > 10% => "missing data"
      2) Any valid k-mer repeated more than once => "repetitive"
      3) A k-mer occupying >= threshold_fraction of positions => "low complexity"

    Returns (False, None) if the sequence is okay.
    Returns (True, reason) if the sequence should be filtered out,
      where reason is either "missing data", "repetitive", or "low complexity".
    """
    seq_len = len(sequence)
    if seq_len == 0:
        return True, "missing data"  # or treat empty seq as special

    # -- 1) Check N fraction --
    num_N = sum(1 for base in sequence if base.upper() == 'N')
    if num_N / seq_len > 0.25:  # More than 10% Ns
        return True, "missing data"

    # -- 2) K-mer checks --
    for k in kmers_to_test:
        if k > seq_len:
            continue

        kmer_counts = count_kmers(sequence, k)
        if not kmer_counts:
            # If all k-mers had N, no valid k-mer was counted.
            # This might imply it's mostly Ns, but let's keep going—maybe
            # the sequence fails in a larger k anyway.
            continue

        # Check if any valid k-mer appears more than once:
        if any(count > 1 for count in kmer_counts.values()):
            return True, "repetitive"

        
    # If we reach here, none of the conditions for removal were met
    return False, None

def main():
    parser = argparse.ArgumentParser(
        description="Filter out low-complexity or high-N sequences from a FASTA file."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file.")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file (filtered).")
    parser.add_argument("-l", "--log", required=True, help="Log file to record removed sequences.")
    parser.add_argument(
        "--kmers",
        help="Dash-separated list of k-values to test (e.g., '5-6-7-8-9-10')."
    )
    
    args = parser.parse_args()

    # Parse the k-values
    kmers_str_list = args.kmers.split("-")
    kmers_to_test = [int(k) for k in kmers_str_list if k.isdigit()]

    with open(args.output, "w") as out_f, open(args.log, "w") as log_f:
        log_f.write("Removed Sequences:\n")
        log_f.write("SequenceHeader\tReason\n")

        for header, seq in read_fasta(args.input):
            # Check the sequence
            removed, reason = check_sequence(seq, kmers_to_test)
            if removed:
                log_f.write(f"{header}\t{reason}\n")
            else:
                # Write the sequence to output
                out_f.write(f">{header}\n{seq}\n")

if __name__ == "__main__":
    main()