# ============================================================
# 05_ortholog_detection.smk - Synteny-Based Ortholog Detection
# ============================================================
#
# Identifies orthologous introner loci across samples using a
# synteny-based approach with flanking region alignment.
#
# Pipeline Overview:
# - Phase 0: GTF to BED conversion for gene annotations
# - Phase 1A: Synteny mapping (upstream/downstream gene context)
# - Phase 1B: Flank extraction and BWA indexing
# - Phase 1C: Initial ortholog pairing (standard alignments)
# - Phase 2: Rescue alignments for Group1↔Group2 pairs
# - Final: Build genotype matrix, fix orientations, annotate
#
# Adapted from: /scratch1/chris/introner-genotyping-pipeline/rules/ortholog_detection_pipeline_v2.smk
#
# ============================================================

import os
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

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

# ============================================================
# PHASE 0: Gene Annotation Processing
# ============================================================

rule gtf_to_gene_bed:
    """
    Convert GTF gene annotations to simplified BED format for bedtools.

    Extracts gene-level features from GTF and outputs 6-column BED format:
    chrom, start (0-based), end (1-based), name (gene_id), score (.), strand
    """
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
    """
    Compute codon-level insertion fingerprints for introner loci.

    For each introner with a gene and splice site annotation, maps the
    insertion site to a (gene_id, codon_number, codon_offset) fingerprint
    using the splice site position and GTF CDS exon structure.

    Used downstream by classify_sharing_status to distinguish ancestral
    shared introners from independent insertions at the same locus.
    """
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
    """
    Detect tandem duplications within each sample's introner candidate loci.

    For each pair of introners in the same sample within a configurable
    genomic window, checks if their body sequences have very high identity.
    Groups them into tandem clusters via union-find.

    Tandem clusters arise from gene duplications and produce introners
    that the flank-based ortholog detection cannot distinguish (because
    their flanking regions are also duplicated). The output is used to
    flag affected ortholog groups in the splitting step.
    """
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


# ============================================================
# PHASE 1A: Synteny-Based Context Mapping
# ============================================================

rule annotate_introner_context:
    """
    Create introner synteny map with three-part genomic fingerprint.

    For each introner locus, identify:
    1. Host gene (containing gene via bedtools intersect)
    2. Upstream gene neighbor (bedtools closest -io -id)
    3. Downstream gene neighbor (bedtools closest -io -iu)

    The synteny fingerprint (Upstream <-- Host --> Downstream) enables
    ortholog detection across samples with different introner content.
    """
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


# ============================================================
# PHASE 1B: Flanking Region Extraction and Indexing
# ============================================================

rule extract_flanks:
    """
    Extract left and right flanking sequences from candidate loci.

    Splits the candidate loci FASTA (which contains introner + flanks)
    into separate left flank and right flank FASTA files for alignment.
    """
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
    """
    Align flanking sequences from query sample to target genome.

    This is the core ortholog detection step: if left and right flanks
    from a query introner align close together in the target genome,
    the locus is likely orthologous.
    """
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


# ============================================================
# PHASE 1C: Initial Ortholog Pairing
# ============================================================

rule pair_orthologs_initial:
    """
    PHASE 1 INITIAL PASS: Pair orthologs from Phase 1 BAM alignments.

    Uses left/right flank alignment positions to identify orthologous
    introner loci across samples. Produces initial results that identify
    scenario 3 (missing data) loci for Phase 2 rescue.

    Scenarios:
    1. Both flanks align together → introner present
    2. Only flanks align (no introner) → introner absent
    3. No/poor alignment → missing data (requires rescue)
    """
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


# ============================================================
# PHASE 2: Rescue Alignments for Distant Pairs
# ============================================================

rule extract_scenario3_flanks:
    """
    Extract flank sequences for scenario 3 loci needing rescue.

    Identifies all scenario 3 (missing data) cases from initial
    ortholog detection and extracts corresponding flank sequences
    for Group1↔Group2 rescue alignments with relaxed parameters.
    """
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
    """
    Align scenario 3 flanks with relaxed BWA parameters.

    Uses relaxed parameters for divergent Group1↔Group2 sequences:
    -k: shorter seed length
    -W: shorter minimum seed length
    -r: less stringent re-seeding
    -A/-B: relaxed scoring (lower mismatch penalty)
    """
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


# ============================================================
# PHASE 1 RESCUE PASS: Incorporate Rescue Data
# ============================================================

rule pair_orthologs_rescue:
    """
    PHASE 1 RESCUE PASS: Pair orthologs with Phase 2 rescue data.

    Re-runs ortholog pairing using BOTH Phase 1 BAM files AND
    Phase 2 rescue BAMs. Fills in scenario 3 loci with rescue data
    to create the final ortholog results.
    """
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


# ============================================================
# FINAL: Build and Annotate Genotype Matrix
# ============================================================

rule build_genotype_matrix:
    """
    Build genotype matrix from ortholog detection results.

    Creates a matrix of introner presence/absence across all samples
    using the network-based approach to group orthologous loci.
    """
    input:
        orthologs = GENOTYPING_DIR / "ortholog_results.tsv"
    output:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.raw.tsv"
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/introner_network.py \
            {input.orthologs} {output.genotype_matrix}
        """


rule fix_orientations:
    """
    Fix introner orientations in the genotype matrix.

    Ensures consistent orientation (relative to gene direction)
    across all samples for each ortholog group.
    """
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


rule classify_sharing_status:
    """
    Classify sharing status for introner ortholog groups.

    Compares codon-level insertion fingerprints within each ortholog group
    to produce two independent classifications:

    within_group_status (consistency within each clade):
      - consistent:  fingerprints agree within each clade
      - discordant:  at least one clade has members at multiple sites
      - uncertain:   insufficient fingerprint data
      - singleton:   only 1 presence=1 member in the group

    cross_group_status (G1 vs G2 ancestry):
      - ancestral:    same insertion site (aa-context match) + same family
      - independent:  different site, different family, or legacy-only match
      - uncertain:    insufficient data for cross-group comparison
      - NA:           no cross-group members
    """
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.oriented.tsv",
        fingerprints = expand(
            GENOTYPING_DIR / "insertion_fingerprints" / "{sample}.fingerprints.tsv",
            sample=ALL_SAMPLES)
    output:
        verified_matrix = GENOTYPING_DIR / "genotype_matrix.verified.tsv",
        summary = GENOTYPING_DIR / "insertion_fingerprints" / "sharing_summary.tsv"
    params:
        codon_tolerance = 3
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/classify_sharing_status.py \
            --matrix {input.genotype_matrix} \
            --fingerprints {input.fingerprints} \
            --output {output.verified_matrix} \
            --summary {output.summary} \
            --codon-tolerance {params.codon_tolerance}
        """


rule split_overmerged_orthologs:
    """
    Split over-merged ortholog groups using codon/intron position clusters.

    The flank-based ortholog detection can incorrectly merge introners
    that are at distinct positions within the same gene because their
    flanking sequences are similar (especially in tandem duplications,
    paralogous gene copies, or genes with internal repeats).

    This rule uses the codon/intron fingerprints to identify ortholog
    groups where Group 1 OR Group 2 contains members at multiple distinct
    positions (within-group over-merge). Such groups are split into
    separate sub-ortholog-groups, one per position cluster.

    Cross-group differences (G1 at one position, G2 at another) are
    NOT split — they represent meaningful biology (independent insertions
    at the same approximate locus).

    The split matrix preserves the original ortholog ID in a new column
    so that the original flank-based grouping can still be referenced.
    """
    input:
        verified_matrix = GENOTYPING_DIR / "genotype_matrix.verified.tsv",
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
            sample=ALL_SAMPLES)
    output:
        split_matrix = GENOTYPING_DIR / "genotype_matrix.split.tsv",
        mapping = GENOTYPING_DIR / "insertion_fingerprints" / "split_mapping.tsv"
    params:
        codon_tolerance = 3
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/split_overmerged_orthologs.py \
            --matrix {input.verified_matrix} \
            --fingerprints {input.fingerprints} \
            --beds {input.beds} \
            --gtfs {input.gtfs} \
            --tandems {input.tandems} \
            --output {output.split_matrix} \
            --mapping {output.mapping} \
            --codon-tolerance {params.codon_tolerance}
        """


rule reclassify_after_split:
    """
    Re-run sharing status classification on the split matrix.

    The splitting changes ortholog group composition, so within-group
    consistency and cross-group classifications must be recomputed.
    """
    input:
        split_matrix = GENOTYPING_DIR / "genotype_matrix.split.tsv",
        fingerprints = expand(
            GENOTYPING_DIR / "insertion_fingerprints" / "{sample}.fingerprints.tsv",
            sample=ALL_SAMPLES)
    output:
        reclassified_matrix = GENOTYPING_DIR / "genotype_matrix.split_verified.tsv",
        summary = GENOTYPING_DIR / "insertion_fingerprints" / "sharing_summary_split.tsv"
    params:
        codon_tolerance = 3
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/classify_sharing_status.py \
            --matrix {input.split_matrix} \
            --fingerprints {input.fingerprints} \
            --output {output.reclassified_matrix} \
            --summary {output.summary} \
            --codon-tolerance {params.codon_tolerance}
        """


rule split_cross_group_mispairs_pass2:
    """
    Second-pass splitting for cross-group mispairs that surfaced only after
    the first reclassify.

    The first pass of `split_overmerged_orthologs` catches within-group
    over-merges plus any cross-group mispairs visible at the initial
    classification. But within-group splits create new sub-groups (_split1,
    _split2) whose cross-group relationship is only determined by the
    subsequent `reclassify_after_split`. If a within-split sub-group now
    contains G1 and G2 members at genuinely different sites, it needs
    another round of cross-group splitting.

    Runs the same splitter script with --cross-group-only, so within-group
    logic is skipped and only cross_group_reason flagged groups get split.
    """
    input:
        split_verified_matrix = GENOTYPING_DIR / "genotype_matrix.split_verified.tsv",
        fingerprints = expand(
            GENOTYPING_DIR / "insertion_fingerprints" / "{sample}.fingerprints.tsv",
            sample=ALL_SAMPLES),
        beds = expand(
            BLAST_DIR / "{sample}.candidate_loci.filtered.bed",
            sample=ALL_SAMPLES),
        gtfs = expand(
            ANNOTATIONS_DIR / "{sample}.gtf",
            sample=ALL_SAMPLES)
    output:
        split_matrix = GENOTYPING_DIR / "genotype_matrix.cross_split.tsv",
        mapping = GENOTYPING_DIR / "insertion_fingerprints" / "cross_split_mapping.tsv"
    params:
        codon_tolerance = 3
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/split_overmerged_orthologs.py \
            --matrix {input.split_verified_matrix} \
            --fingerprints {input.fingerprints} \
            --beds {input.beds} \
            --gtfs {input.gtfs} \
            --output {output.split_matrix} \
            --mapping {output.mapping} \
            --codon-tolerance {params.codon_tolerance} \
            --cross-group-only
        """


rule reclassify_after_cross_split:
    """
    Re-run sharing-status classification on the post-pass-2 matrix so the
    new cgsplit sub-groups get `within_group_status = consistent` and
    `cross_group_status = NA` labels reflecting their clade-specific
    composition.
    """
    input:
        split_matrix = GENOTYPING_DIR / "genotype_matrix.cross_split.tsv",
        fingerprints = expand(
            GENOTYPING_DIR / "insertion_fingerprints" / "{sample}.fingerprints.tsv",
            sample=ALL_SAMPLES)
    output:
        reclassified_matrix = GENOTYPING_DIR / "genotype_matrix.cross_split_verified.tsv",
        summary = GENOTYPING_DIR / "insertion_fingerprints" / "sharing_summary_cross_split.tsv"
    params:
        codon_tolerance = 3
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/classify_sharing_status.py \
            --matrix {input.split_matrix} \
            --fingerprints {input.fingerprints} \
            --output {output.reclassified_matrix} \
            --summary {output.summary} \
            --codon-tolerance {params.codon_tolerance}
        """


rule compare_introner_sequences:
    """
    Refine within_group_status and cross_group_status using body sequence identity.

    Extracts introner body sequences from indexed genome FASTAs and computes
    pairwise identity. Overwrites the codon-level status columns with refined
    values using two thresholds:
      - within-group identity (default 0.80): resolves uncertain intergenic
        groups to consistent (if same family + high identity), flags low_identity
      - cross-group identity (default 0.60): resolves uncertain cross-group
        cases to likely_ancestral/likely_independent

    Adds audit columns: within_group_identity, cross_group_identity.
    """
    input:
        verified_matrix = GENOTYPING_DIR / "genotype_matrix.cross_split_verified.tsv",
        genome_indices = expand(
            ASSEMBLIES_DIR / "{sample}.vg_paths.fa.fai",
            sample=ALL_SAMPLES)
    output:
        seq_verified_matrix = GENOTYPING_DIR / "genotype_matrix.seq_verified.tsv",
        summary = GENOTYPING_DIR / "insertion_fingerprints" / "sequence_comparison_summary.tsv"
    params:
        genome_dir = ASSEMBLIES_DIR,
        flanking_length = FLANK_LENGTH,
        within_threshold = 0.80,
        cross_threshold = 0.60
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/compare_introner_sequences.py \
            --matrix {input.verified_matrix} \
            --genome-dir {params.genome_dir} \
            --output {output.seq_verified_matrix} \
            --summary {output.summary} \
            --flanking-length {params.flanking_length} \
            --within-group-identity-threshold {params.within_threshold} \
            --cross-group-identity-threshold {params.cross_threshold}
        """


rule annotate_missing_data:
    """
    Annotate missing gene and family data in the genotype matrix.

    Strategy:
    1. Drop ortholog groups with within_group_status of 'discordant' or
       'uncertain' (unusable for pi/dxy calculations)

    2. Gene annotation (ALL rows): Use bedtools overlap to find genes
       at each introner's genomic coordinates (works for all scenarios)

    3. Family annotation (Scenario 1 only): Use sequence similarity
       to match against reference introner families (presence=1 only)

    Output: Fully annotated genotype matrix ready for downstream analysis
    """
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.seq_verified.tsv",
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


# ============================================================
# TARGET RULES
# ============================================================

rule all_ortholog_detection:
    """
    Target: Complete ortholog detection pipeline
    """
    input:
        GENOTYPING_DIR / "genotype_matrix.tsv"


rule ortholog_alignments_only:
    """
    Target: Generate all flank alignments without building matrix
    """
    input:
        [ORTHOLOG_DIR / f"{query}_{target}.left_flank.sorted.bam"
         for query in ALL_SAMPLES
         for target in ALL_SAMPLES if target != query],
        [ORTHOLOG_DIR / f"{query}_{target}.right_flank.sorted.bam"
         for query in ALL_SAMPLES
         for target in ALL_SAMPLES if target != query]
