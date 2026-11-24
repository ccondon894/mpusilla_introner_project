import os
from pathlib import Path

# Output directories
ORTHOLOG_DIR = GENOTYPING_DIR / "ortholog_detection"
ORTHOLOG_LOG_DIR = ORTHOLOG_DIR / "logs"
PROCESSED_ANN_DIR = ORTHOLOG_DIR / "processed_annotations"
RESCUE_DIR = ORTHOLOG_DIR / "rescue_alignments"

# BLAST results directory (input from 04_blast_genotyping.smk)
BLAST_DIR = GENOTYPING_DIR / "blast_results"

# Reference data for family annotation
REF_INTRONER_DIR = Path(config["paths"].get("ref_introner_data",
    "/scratch1/chris/introner-genotyping-pipeline/introner_files_for_Github/"))

# Generate list of distant pairs for Phase 2 rescue
# Group1 ↔ Group2 alignments need relaxed parameters
DISTANT_PAIRS = (
    [(q, t) for q in GROUP1_SAMPLES for t in GROUP2_SAMPLES] +
    [(q, t) for q in GROUP2_SAMPLES for t in GROUP1_SAMPLES]
)

# Helper function for cross-sample alignments
def get_targets(query):
    """Return all samples except the query"""
    return [s for s in ALL_SAMPLES if s != query]

# BWA parameters from config
BWA_THREADS = config["params"]["bwa"]["threads"]
BWA_RESCUE_K = config["params"]["bwa_rescue"]["k"]
BWA_RESCUE_W = config["params"]["bwa_rescue"]["W"]
BWA_RESCUE_R = config["params"]["bwa_rescue"]["r"]
BWA_RESCUE_A = config["params"]["bwa_rescue"]["A"]
BWA_RESCUE_B = config["params"]["bwa_rescue"]["B"]

# Annotation parameters
FLANK_LENGTH = config["params"]["flanks"]["extraction_length"]
SIMILARITY_CUTOFF = config["params"]["similarity"]["cutoff"]

rule gtf_to_gene_bed:
    """Convert GTF gene annotations to simplified BED format for bedtools."""
    input:
        gtf = ANNOTATIONS_DIR / "{sample}.gtf"
    output:
        bed = PROCESSED_ANN_DIR / "{sample}.gene.bed"
    shell:
        """
        mkdir -p {PROCESSED_ANN_DIR}
        awk 'BEGIN{{OFS="\\t"}} $3 == "gene" {{
            gene_id = "";
            for (i=9; i<=NF; i++) {{
                if ($i == "gene_id") {{
                    if (i+1 <= NF) {{
                        val = $(i+1);
                        gsub(/[";]/, "", val);
                        gene_id = val;
                    }}
                    break;
                }}
            }}
            print $1, $4-1, $5, gene_id, ".", $7
        }}' {input.gtf} > {output.bed}
        """


rule compute_insertion_fingerprints:
    """Compute codon-level insertion fingerprints for introner loci."""
    input:
        bed = BLAST_DIR / "{sample}.candidate_loci.filtered.bed",
        gtf = ANNOTATIONS_DIR / "{sample}.gtf",
        genome = ASSEMBLIES_DIR / "{sample}.vg_paths.fa",
        genome_index = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.fai"
    output:
        fingerprints = GENOTYPING_DIR / "insertion_fingerprints" / "{sample}.fingerprints.tsv"
    shell:
        """
        mkdir -p {GENOTYPING_DIR}/insertion_fingerprints
        python {PROJECT_ROOT}/scripts/genotyping/compute_insertion_fingerprints.py \
            --bed {input.bed} \
            --gtf {input.gtf} \
            --genome {input.genome} \
            --sample {wildcards.sample} \
            --output {output.fingerprints}
        """


rule detect_tandem_duplicates:
    """Detect tandem duplications within each sample's introner candidate loci."""
    input:
        bed = BLAST_DIR / "{sample}.candidate_loci.filtered.bed",
        genome = ASSEMBLIES_DIR / "{sample}.vg_paths.fa",
        genome_index = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.fai"
    output:
        tandems = GENOTYPING_DIR / "insertion_fingerprints" / "{sample}.tandems.tsv"
    params:
        window_bp = 10000,
        min_identity = 0.95
    shell:
        """
        mkdir -p {GENOTYPING_DIR}/insertion_fingerprints
        python {PROJECT_ROOT}/scripts/genotyping/detect_tandem_duplicates.py \
            --bed {input.bed} \
            --genome {input.genome} \
            --sample {wildcards.sample} \
            --output {output.tandems} \
            --window-bp {params.window_bp} \
            --min-identity {params.min_identity}
        """

# PHASE 1A: Synteny-Based Context Mapping
rule annotate_introner_context:
    """Create introner synteny map with three-part genomic fingerprint."""
    input:
        introner_bed = BLAST_DIR / "{sample}.candidate_loci.filtered.bed",
        gene_bed = PROCESSED_ANN_DIR / "{sample}.gene.bed"
    output:
        intersect_genes = ORTHOLOG_DIR / "{sample}.intersect_genes.txt",
        upstream_neighbors = ORTHOLOG_DIR / "{sample}.upstream_neighbors.txt",
        downstream_neighbors = ORTHOLOG_DIR / "{sample}.downstream_neighbors.txt",
        synteny_map = ORTHOLOG_DIR / "{sample}.introner_synteny_map.tsv"
    shell:
        """
        mkdir -p {ORTHOLOG_DIR}

        # Sort BED files for bedtools
        sorted_genes=$(mktemp)
        sorted_introner=$(mktemp)

        sort -k1,1 -k2,2n {input.gene_bed} > $sorted_genes
        sort -k1,1 -k2,2n {input.introner_bed} > $sorted_introner

        # Step 1: Find genes containing/overlapping introners
        bedtools intersect \
            -a $sorted_introner \
            -b $sorted_genes \
            -wo > {output.intersect_genes}

        # Step 2: Find UPSTREAM neighbor (ignoring self-hits)
        bedtools closest \
            -a $sorted_genes \
            -b $sorted_genes \
            -io \
            -id \
            -D ref > {output.upstream_neighbors}

        # Step 3: Find DOWNSTREAM neighbor (ignoring self-hits)
        bedtools closest \
            -a $sorted_genes \
            -b $sorted_genes \
            -io \
            -iu \
            -D ref > {output.downstream_neighbors}

        # Step 4: Combine into synteny map
        python {PROJECT_ROOT}/scripts/genotyping/create_introner_context.py \
            {output.intersect_genes} \
            {output.upstream_neighbors} \
            {output.downstream_neighbors} \
            {output.synteny_map}

        # Clean up temp files
        rm $sorted_genes $sorted_introner
        """


# PHASE 1B: Flanking Region Extraction and Indexing
rule extract_flanks:
    """Extract left and right flanking sequences from candidate loci."""
    input:
        fa = BLAST_DIR / "{sample}.candidate_loci.filtered.fa"
    output:
        left = ORTHOLOG_DIR / "{sample}.left_flanks.fa",
        right = ORTHOLOG_DIR / "{sample}.right_flanks.fa"
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/extract_flanks.py \
            {input.fa} {output.left} {output.right}
        """


rule bwa_align_flanks:
    """Align flanking sequences from query sample to target genome."""
    input:
        query_left = ORTHOLOG_DIR / "{query}.left_flanks.fa",
        query_right = ORTHOLOG_DIR / "{query}.right_flanks.fa",
        target_genome = ASSEMBLIES_DIR / "{target}.vg_paths.fa",
        target_index = ASSEMBLIES_DIR / "{target}.vg_paths.fa.bwt"
    output:
        left_bam = ORTHOLOG_DIR / "{query}_{target}.left_flank.sorted.bam",
        right_bam = ORTHOLOG_DIR / "{query}_{target}.right_flank.sorted.bam",
        left_bai = ORTHOLOG_DIR / "{query}_{target}.left_flank.sorted.bam.bai",
        right_bai = ORTHOLOG_DIR / "{query}_{target}.right_flank.sorted.bam.bai"
    threads: BWA_THREADS
    shell:
        """
        # Align left flanks
        bwa mem -t {threads} {input.target_genome} {input.query_left} | \
        samtools sort -@ {threads} -o {output.left_bam}
        samtools index {output.left_bam}

        # Align right flanks
        bwa mem -t {threads} {input.target_genome} {input.query_right} | \
        samtools sort -@ {threads} -o {output.right_bam}
        samtools index {output.right_bam}
        """

# PHASE 1C: Initial Ortholog Pairing
rule pair_orthologs_initial:
    """PHASE 1 INITIAL PASS: Pair orthologs from Phase 1 BAM alignments."""
    input:
        left_bam = [ORTHOLOG_DIR / f"{query}_{target}.left_flank.sorted.bam"
                    for query in ALL_SAMPLES
                    for target in ALL_SAMPLES if target != query],
        right_bam = [ORTHOLOG_DIR / f"{query}_{target}.right_flank.sorted.bam"
                     for query in ALL_SAMPLES
                     for target in ALL_SAMPLES if target != query],
        bed_files = [BLAST_DIR / f"{sample}.candidate_loci.filtered.bed"
                     for sample in ALL_SAMPLES]
    output:
        orthologs = GENOTYPING_DIR / "ortholog_results.initial.tsv"
    params:
        ortholog_dir = ORTHOLOG_DIR,
        blast_dir = BLAST_DIR,
        assembly_dir = ASSEMBLIES_DIR,
        processed_ann_dir = PROCESSED_ANN_DIR
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/introner_genotype_matrix_builder.py \
            {params.ortholog_dir} {output.orthologs} --mode initial \
            --bed-dir {params.blast_dir} \
            --genome-dir {params.assembly_dir} \
            --processed-ann-dir {params.processed_ann_dir}
        """


# PHASE 2: Rescue Alignments for Distant Pairs

rule extract_scenario3_flanks:
    """Extract flank sequences for scenario 3 loci needing rescue."""
    input:
        ortholog_results = GENOTYPING_DIR / "ortholog_results.initial.tsv",
        left_flanks = [ORTHOLOG_DIR / f"{sample}.left_flanks.fa" for sample in ALL_SAMPLES],
        right_flanks = [ORTHOLOG_DIR / f"{sample}.right_flanks.fa" for sample in ALL_SAMPLES]
    output:
        rescue_flanks = [RESCUE_DIR / f"{query}_to_{target}.scenario3_flanks.fa"
                         for query, target in DISTANT_PAIRS]
    params:
        ortholog_dir = ORTHOLOG_DIR,
        rescue_dir = RESCUE_DIR
    shell:
        """
        mkdir -p {params.rescue_dir}
        python {PROJECT_ROOT}/scripts/genotyping/extract_scenario3_flanks.py \
            {input.ortholog_results} {params.ortholog_dir} {params.rescue_dir}
        """


rule align_rescue_flanks:
    """Align scenario 3 flanks with relaxed BWA parameters."""
    input:
        rescue_flanks = RESCUE_DIR / "{query}_to_{target}.scenario3_flanks.fa",
        target_genome = ASSEMBLIES_DIR / "{target}.vg_paths.fa",
        target_index = ASSEMBLIES_DIR / "{target}.vg_paths.fa.bwt"
    output:
        bam = RESCUE_DIR / "{query}_to_{target}.rescue.sorted.bam",
        bai = RESCUE_DIR / "{query}_to_{target}.rescue.sorted.bam.bai"
    params:
        k = BWA_RESCUE_K,
        W = BWA_RESCUE_W,
        r = BWA_RESCUE_R,
        A = BWA_RESCUE_A,
        B = BWA_RESCUE_B
    threads: BWA_THREADS
    shell:
        """
        # Relaxed BWA parameters for divergent sequences
        bwa mem -t {threads} \
            -k {params.k} -W {params.W} -r {params.r} \
            -A {params.A} -B {params.B} \
            {input.target_genome} {input.rescue_flanks} | \
        samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam}
        """


# PHASE 1 RESCUE PASS: Incorporate Rescue Data

rule pair_orthologs_rescue:
    """PHASE 1 RESCUE PASS: Pair orthologs with Phase 2 rescue data."""
    input:
        left_bam = [ORTHOLOG_DIR / f"{query}_{target}.left_flank.sorted.bam"
                    for query in ALL_SAMPLES
                    for target in ALL_SAMPLES if target != query],
        right_bam = [ORTHOLOG_DIR / f"{query}_{target}.right_flank.sorted.bam"
                     for query in ALL_SAMPLES
                     for target in ALL_SAMPLES if target != query],
        rescue_bams = [RESCUE_DIR / f"{query}_to_{target}.rescue.sorted.bam"
                       for query, target in DISTANT_PAIRS],
        bed_files = [BLAST_DIR / f"{sample}.candidate_loci.filtered.bed"
                     for sample in ALL_SAMPLES]
    output:
        orthologs = GENOTYPING_DIR / "ortholog_results.tsv"
    params:
        ortholog_dir = ORTHOLOG_DIR,
        blast_dir = BLAST_DIR,
        assembly_dir = ASSEMBLIES_DIR,
        processed_ann_dir = PROCESSED_ANN_DIR
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/introner_genotype_matrix_builder.py \
            {params.ortholog_dir} {output.orthologs} --mode rescue \
            --bed-dir {params.blast_dir} \
            --genome-dir {params.assembly_dir} \
            --processed-ann-dir {params.processed_ann_dir}
        """


# FINAL: Build and Annotate Genotype Matrix
rule build_genotype_matrix:
    """Build genotype matrix from ortholog detection results."""
    input:
        orthologs = GENOTYPING_DIR / "ortholog_results.tsv"
    output:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.raw.tsv",
        locus_group_members = GENOTYPING_DIR / "locus_group_members.tsv"
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/introner_network.py \
            {input.orthologs} {output.genotype_matrix} {output.locus_group_members}
        """


rule fix_orientations:
    """Fix introner orientations in the genotype matrix."""
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.raw.tsv"
    output:
        fixed_matrix = GENOTYPING_DIR / "genotype_matrix.oriented.tsv"
    params:
        assembly_dir = ASSEMBLIES_DIR
    threads: config["runtime"]["threads_max"]
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/check_orientation.py \
            {input.genotype_matrix} {params.assembly_dir} {output.fixed_matrix} \
            --threads {threads}
        """


rule resolve_ortholog_groups:
    """Consolidated ortholog resolver."""
    input:
        resolver_script = PROJECT_ROOT / "scripts" / "genotyping" / "resolve_ortholog_groups.py",
        classifier_script = PROJECT_ROOT / "scripts" / "genotyping" / "classify_sharing_status.py",
        splitter_script = PROJECT_ROOT / "scripts" / "genotyping" / "split_overmerged_orthologs.py",
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.oriented.tsv",
        locus_group_members = GENOTYPING_DIR / "locus_group_members.tsv",
        fingerprints = expand(
            GENOTYPING_DIR / "insertion_fingerprints" / "{sample}.fingerprints.tsv",
            sample=ALL_SAMPLES),
        beds = expand(
            BLAST_DIR / "{sample}.candidate_loci.filtered.bed",
            sample=ALL_SAMPLES),
        gtfs = expand(
            ANNOTATIONS_DIR / "{sample}.gtf",
            sample=ALL_SAMPLES),
        tandems = expand(
            GENOTYPING_DIR / "insertion_fingerprints" / "{sample}.tandems.tsv",
            sample=ALL_SAMPLES),
        genome_indices = expand(
            ASSEMBLIES_DIR / "{sample}.vg_paths.fa.fai",
            sample=ALL_SAMPLES)
    output:
        resolved_matrix = GENOTYPING_DIR / "genotype_matrix.resolved.tsv",
        insertion_group_map = GENOTYPING_DIR / "insertion_fingerprints" / "insertion_group_map.tsv",
        independent_insertion_events = GENOTYPING_DIR / "insertion_fingerprints" / "independent_insertion_events.tsv",
        resolution_summary = GENOTYPING_DIR / "insertion_fingerprints" / "ortholog_resolution_summary.tsv"
    params:
        genome_dir = ASSEMBLIES_DIR,
        codon_tolerance = 3,
        flanking_length = FLANK_LENGTH,
        within_threshold = 0.80,
        cross_threshold = 0.60
    shell:
        """
        python {input.resolver_script} \
            --matrix {input.genotype_matrix} \
            --locus-group-members {input.locus_group_members} \
            --fingerprints {input.fingerprints} \
            --beds {input.beds} \
            --gtfs {input.gtfs} \
            --tandems {input.tandems} \
            --genome-dir {params.genome_dir} \
            --output {output.resolved_matrix} \
            --insertion-group-map {output.insertion_group_map} \
            --independent-insertion-events {output.independent_insertion_events} \
            --resolution-summary {output.resolution_summary} \
            --codon-tolerance {params.codon_tolerance} \
            --flanking-length {params.flanking_length} \
            --within-group-identity-threshold {params.within_threshold} \
            --cross-group-identity-threshold {params.cross_threshold}
        """


rule annotate_missing_data:
    """Annotate missing gene and family data in the genotype matrix."""
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.resolved.tsv",
        gene_beds = [PROCESSED_ANN_DIR / f"{sample}.gene.bed" for sample in ALL_SAMPLES],
        fasta_files = [BLAST_DIR / f"{sample}.candidate_loci.filtered.fa" for sample in ALL_SAMPLES]
    output:
        annotated_matrix = GENOTYPING_DIR / "genotype_matrix.tsv",
        log_file = ORTHOLOG_LOG_DIR / "genotype_matrix_annotation.log"
    params:
        processed_ann_dir = PROCESSED_ANN_DIR,
        fasta_dir = BLAST_DIR,
        assembly_dir = ASSEMBLIES_DIR,
        ref_dir = REF_INTRONER_DIR,
        flank_length = FLANK_LENGTH,
        similarity_cutoff = SIMILARITY_CUTOFF,
        annotation_method = "ortholog_group"
    shell:
        """
        mkdir -p {ORTHOLOG_LOG_DIR}
        python {PROJECT_ROOT}/scripts/genotyping/annotate_missing_introners.py \
            --genotype_matrix {input.genotype_matrix} \
            --processed_ann_dir {params.processed_ann_dir} \
            --fasta_dir {params.fasta_dir} \
            --genome_dir {params.assembly_dir} \
            --ref_dir {params.ref_dir} \
            --output {output.annotated_matrix} \
            --output_log {output.log_file} \
            --flanking_length {params.flank_length} \
            --similarity_cutoff {params.similarity_cutoff} \
            --annotation_method {params.annotation_method}
        """


rule all_ortholog_detection:
    """Target: Complete ortholog detection pipeline"""
    input:
        GENOTYPING_DIR / "genotype_matrix.tsv"


rule ortholog_alignments_only:
    """Target: Generate all flank alignments without building matrix"""
    input:
        [ORTHOLOG_DIR / f"{query}_{target}.left_flank.sorted.bam"
         for query in ALL_SAMPLES
         for target in ALL_SAMPLES if target != query],
        [ORTHOLOG_DIR / f"{query}_{target}.right_flank.sorted.bam"
         for query in ALL_SAMPLES
         for target in ALL_SAMPLES if target != query]
