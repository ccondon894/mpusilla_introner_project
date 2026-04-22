# ============================================================
# 04_blast_genotyping.smk - BLAST-based Introner Detection
# ============================================================
#
# Identifies candidate introner loci by BLAST homology search
# against the introner truth set.
#
# Steps:
# 1. Create BLAST database for each sample genome
# 2. BLAST introner truth set against sample genomes
# 3. Filter hits by identity and length coverage
# 4. Convert to BED format and filter
# 5. Extract candidate loci sequences with flanking regions
# 6. Filter low-complexity sequences
# 7. Validate candidates against reference annotations
#
# Adapted from: /scratch1/chris/introner-genotyping-pipeline/rules/blast_genotyping_graph.smk
#
# ============================================================

import os
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Input files
INTRONER_TRUTH_SET = PROJECT_ROOT / "resources" / "introner_truth_set.fa"

# Output directories
BLAST_DIR = GENOTYPING_DIR / "blast_results"
BLAST_LOG_DIR = BLAST_DIR / "logs"

# Reference data for validation
REF_DIR = Path(config["paths"].get("ref_introner_data",
    "/scratch1/chris/introner-genotyping-pipeline/introner_files_for_Github/"))

# Structural analysis configuration (disabled by default)
ENABLE_STRUCTURAL_ANALYSIS = False
VALIDATE_MICROMONAS_PATTERNS = False

# Parameters from config
BLAST_IDENTITY = config["params"]["blast"]["identity_threshold"]
BLAST_THREADS = config["params"]["blast"]["threads"]
SIMILARITY_CUTOFF = config["params"]["similarity"]["cutoff"]
FLANK_LENGTH = config["params"]["flanks"]["extraction_length"]

# ============================================================
# RULES
# ============================================================

rule blast_introners:
    """
    BLAST introner truth set against sample genome
    """
    input:
        db = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.nhr",
        sequences = INTRONER_TRUTH_SET
    output:
        results = BLAST_DIR / "{sample}_blast_results.txt"
    params:
        db_prefix = lambda wildcards: str(ASSEMBLIES_DIR / f"{wildcards.sample}.vg_paths.fa")
    threads: BLAST_THREADS
    shell:
        """
        mkdir -p {BLAST_DIR}
        blastn -num_threads {threads} \
            -query {input.sequences} \
            -db {params.db_prefix} \
            -out {output.results} \
            -outfmt "6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue"
        """


rule filter_blast_results:
    """
    Filter BLAST results by identity (>=80%) and length coverage (80-120%)
    """
    input:
        blast_results = BLAST_DIR / "{sample}_blast_results.txt"
    output:
        filtered = BLAST_DIR / "{sample}_filtered_blast_results.txt"
    params:
        identity = BLAST_IDENTITY,
        min_cov = config["params"]["blast"]["length_coverage_min"],
        max_cov = config["params"]["blast"]["length_coverage_max"]
    shell:
        """
        awk '($3 >= {params.identity}) && ($6 >= {params.min_cov} * $4) && ($6 <= {params.max_cov} * $4)' \
            {input.blast_results} > {output.filtered}
        """


rule convert_blast_to_bed:
    """
    Convert filtered BLAST results to BED format
    """
    input:
        blast = BLAST_DIR / "{sample}_filtered_blast_results.txt"
    output:
        bed = BLAST_DIR / "{sample}_filtered_blast_results.bed"
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} {{
            start = ($11 < $12) ? $11 : $12
            end = ($11 > $12) ? $11 : $12
            print $2, start-1, end, $1, $3, $4, $5, $6, $7, $8, $9, $10, $13
        }}' {input.blast} > {output.bed}
        """


rule assess_blast_results:
    """
    Assess and deduplicate BLAST hits
    """
    input:
        bed = BLAST_DIR / "{sample}_filtered_blast_results.bed"
    output:
        bed = BLAST_DIR / "{sample}.candidate_loci.bed"
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/assess_blast_results_faster.py \
            {input.bed} {output.bed}
        """


rule add_flanks_to_bed:
    """
    Add flanking regions (100bp each side) to candidate loci
    """
    input:
        bed = BLAST_DIR / "{sample}.candidate_loci.bed"
    output:
        bed = BLAST_DIR / "{sample}.candidate_loci_plus_flanks.bed"
    params:
        flank = FLANK_LENGTH
    shell:
        """
        awk '{{print $1, $2 - {params.flank}, $3 + {params.flank}, $4}}' {input.bed} > {output.bed}
        """


rule extract_candidate_sequences:
    """
    Extract FASTA sequences for candidate loci with flanks
    """
    input:
        bed = BLAST_DIR / "{sample}.candidate_loci_plus_flanks.bed",
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        fa = BLAST_DIR / "{sample}.candidate_loci_plus_flanks.fa"
    shell:
        """
        awk 'OFS="\\t" {{print $1, $2, $3, $1":"$2"-"$3"|"$4}}' {input.bed} \
        | bedtools getfasta -fi {input.fa} -bed - -name \
        | sed '/^>/ s/::.*$//' > {output.fa}
        """


rule filter_low_complexity:
    """
    Remove low-complexity sequences from candidates
    """
    input:
        fa = BLAST_DIR / "{sample}.candidate_loci_plus_flanks.fa"
    output:
        fa = BLAST_DIR / "{sample}.candidate_loci_plus_flanks.filtered.fa"
    log:
        BLAST_LOG_DIR / "{sample}.removed_candidate_loci.log"
    shell:
        """
        mkdir -p {BLAST_LOG_DIR}
        python {PROJECT_ROOT}/scripts/genotyping/remove_low_complexity_sequences.py \
            --input {input.fa} \
            --output {output.fa} \
            --log {log} \
            --kmers 25
        """


rule validate_introner_candidates:
    """
    Validate introner candidates using BLAST-based family classification.

    Scores each candidate against all individual reference introner sequences
    (not a single consensus per family), providing more reliable family
    assignments and explicit confidence margins.
    """
    input:
        fa = BLAST_DIR / "{sample}.candidate_loci_plus_flanks.filtered.fa",
        gtf = ANNOTATIONS_DIR / "{sample}.gtf"
    output:
        fa = BLAST_DIR / "{sample}.candidate_loci.validated.fa",
        log1 = BLAST_LOG_DIR / "{sample}.introner_similarity_check.log",
        log2 = BLAST_LOG_DIR / "{sample}.introner_similarity_check_summary.log"
    params:
        ref_dir = REF_DIR,
        similarity = SIMILARITY_CUTOFF
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/validate_introners_blast.py \
            --input_fasta {input.fa} \
            --gtf_file {input.gtf} \
            --ref_dir {params.ref_dir} \
            --output_fasta {output.fa} \
            --output_log {output.log1} \
            --identity_cutoff {params.similarity} \
            --coverage_cutoff 0.8
        """


rule filter_boundary_spanning_introners:
    """
    Filter introner candidates by gene-boundary overlap.

    Rejects candidates whose body spans multiple genes or extends past a
    single gene's transcript boundary (with a small tolerance for extraction
    noise). A biologically valid introner should be either fully inside a
    gene (intronic/exonic/UTR) or fully in an intergenic region.

    Candidates that straddle gene boundaries are most likely false positive
    BLAST hits where a sequence with moderate similarity to an introner
    reference happens to align across a gene boundary.

    Each kept candidate gets a 'genomic_context' tag (intra_gene or
    intergenic) added to its FASTA header. Rejected candidates are logged
    to a separate file for auditing.
    """
    input:
        fa = BLAST_DIR / "{sample}.candidate_loci.validated.fa",
        gtf = ANNOTATIONS_DIR / "{sample}.gtf"
    output:
        fa = BLAST_DIR / "{sample}.candidate_loci.filtered.fa",
        rejections = BLAST_LOG_DIR / "{sample}.boundary_rejections.tsv"
    params:
        flank_length = 100,
        tolerance = 5
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/filter_boundary_spanning_introners.py \
            --input {input.fa} \
            --gtf {input.gtf} \
            --output {output.fa} \
            --log {output.rejections} \
            --flanking-length {params.flank_length} \
            --tolerance {params.tolerance}
        """


rule create_candidate_bed:
    """
    Create BED file from the boundary-filtered candidate FASTA.

    Includes a genomic_context column indicating whether the candidate is
    inside a gene (intra_gene) or in an intergenic region (intergenic).
    """
    input:
        fa = BLAST_DIR / "{sample}.candidate_loci.filtered.fa"
    output:
        bed = BLAST_DIR / "{sample}.candidate_loci.filtered.bed"
    run:
        with open(input.fa, 'r') as f, open(output.bed, 'w') as o:
            o.write("#chrom\tstart\tend\tname\tfamily\tsimilarity\torientation\t"
                    "gene_info\tsplice_info\tgenomic_context\n")

            for line in f:
                if line.startswith(">"):
                    try:
                        line = line[1:].strip()
                        parts = line.split(" ")
                        region_and_id = parts[0]
                        metadata = parts[1:] if len(parts) > 1 else []

                        # Parse region and identifier
                        region, identifier = region_and_id.split("|")
                        chrom, coords = region.split(":")
                        start, end = coords.split("-")

                        # Parse metadata
                        metadata_dict = {}
                        for item in metadata:
                            if "=" in item:
                                key, value = item.split("=", 1)
                                metadata_dict[key] = value

                        family = metadata_dict.get('family', 'unknown')
                        similarity = metadata_dict.get('similarity', 'NA')
                        orientation = metadata_dict.get('orientation', 'NA')
                        gene_info = metadata_dict.get('gene', 'NA')
                        splice_info = metadata_dict.get('splice_sites', metadata_dict.get('splice_site', 'NA'))
                        genomic_context = metadata_dict.get('genomic_context', 'NA')

                        o.write(f"{chrom}\t{start}\t{end}\t{identifier}\t{family}\t{similarity}\t"
                                f"{orientation}\t{gene_info}\t{splice_info}\t{genomic_context}\n")

                    except (ValueError, IndexError) as e:
                        print(f"Warning: Error processing FASTA header: {line}")
                        continue


rule all_blast_genotyping:
    """
    Target: Generate all BLAST candidate loci
    """
    input:
        expand(BLAST_DIR / "{sample}.candidate_loci.filtered.fa", sample=ALL_SAMPLES),
        expand(BLAST_DIR / "{sample}.candidate_loci.filtered.bed", sample=ALL_SAMPLES)
