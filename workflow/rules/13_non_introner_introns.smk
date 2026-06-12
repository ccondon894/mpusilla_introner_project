# ============================================================================
# 13_non_introner_introns.smk
# Non-introner intron presence/absence polymorphism and allele frequency spectrum
#
# Tests whether non-introner introns show presence/absence polymorphism,
# and computes folded + unfolded AFS for comparison with introner dynamics.
# Uses minimap2 alignment for orthologous intron matching across assemblies,
# with coverage validation to eliminate assembly-derived artifacts.
#
# v3: Adds short-read coverage validation via introner-caller --format simple
# ============================================================================

NON_INTRONER_DIR = EVOLUTION_DIR / "non_introner_introns"
NON_INTRONER_COV_DIR = NON_INTRONER_DIR / "coverage"
NON_INTRONER_CIRCOS_DIR = NON_INTRONER_DIR / "circos"
NON_INTRONER_CIRCOS_COV_DIR = NON_INTRONER_CIRCOS_DIR / "coverage"
NON_INTRONER_SCRIPTS = str(PROJECT_ROOT / "scripts" / "evolution")

# Non-reference samples for coverage analysis (CCMP1545 stays PRESENT by definition)
NON_REF_SAMPLES = [s for s in GROUP1_SAMPLES + GROUP2_SAMPLES if s != REFERENCE]

# Reuse existing BAMs and coverage parameters from 06_coverage_calling.smk
# BAM_DIR = GENOTYPING_DIR / "coverage" / "bams"  (defined in 06_coverage_calling.smk)
# MAPQ_THRESHOLD (defined in 06_coverage_calling.smk)


rule build_reference_intron_catalog:
    """
    Build catalog of reference non-introner introns.
    Extracts introns from reference GTF, removes introner-overlapping ones.
    """
    input:
        reference_gtf = ANNOTATIONS_DIR / f"{REFERENCE}.gtf",
        introner_loci = GENOTYPING_DIR / "introner_loci.bed",
    output:
        catalog = NON_INTRONER_DIR / "reference_intron_catalog.tsv",
    conda: "../envs/non_introner_introns.yaml"
    shell:
        """
        mkdir -p {NON_INTRONER_DIR} && \
        python {NON_INTRONER_SCRIPTS}/build_reference_intron_catalog.py \
            --reference_gtf {input.reference_gtf} \
            --introner_loci {input.introner_loci} \
            --output {output.catalog}
        """


rule build_reference_intron_catalog_with_mating:
    """
    Build a Circos-only catalog of reference non-introner introns.
    This keeps the scaffold 2 mating-type region so chromosome density tracks
    show the full intron landscape.
    """
    input:
        reference_gtf = ANNOTATIONS_DIR / f"{REFERENCE}.gtf",
        introner_loci = GENOTYPING_DIR / "introner_loci.bed",
    output:
        catalog = NON_INTRONER_CIRCOS_DIR / "reference_intron_catalog.with_mating.tsv",
    conda: "../envs/non_introner_introns.yaml"
    shell:
        """
        mkdir -p {NON_INTRONER_CIRCOS_DIR} && \
        python {NON_INTRONER_SCRIPTS}/build_reference_intron_catalog.py \
            --reference_gtf {input.reference_gtf} \
            --introner_loci {input.introner_loci} \
            --include_mating_region \
            --output {output.catalog}
        """


rule make_intron_coverage_beds_with_mating:
    """
    Create Circos-only per-sample BED files from the with-mating intron catalog.
    """
    input:
        gtf = lambda wildcards: get_gtf(wildcards.sample),
        catalog = NON_INTRONER_CIRCOS_DIR / "reference_intron_catalog.with_mating.tsv",
        assembly = ASSEMBLIES_DIR / "{sample}.vg_paths.fa",
    output:
        bed = NON_INTRONER_CIRCOS_COV_DIR / "{sample}.intron_loci.bed",
    params:
        sample = "{sample}",
    conda: "../envs/non_introner_introns.yaml"
    shell:
        """
        mkdir -p {NON_INTRONER_CIRCOS_COV_DIR} && \
        python {NON_INTRONER_SCRIPTS}/make_intron_coverage_beds.py \
            --gtf {input.gtf} \
            --reference_catalog {input.catalog} \
            --assembly {input.assembly} \
            --sample {params.sample} \
            --output {output.bed}
        """


rule calculate_intron_locus_depth_with_mating:
    """
    Calculate read depth at Circos-only with-mating intron loci.
    """
    input:
        bed = NON_INTRONER_CIRCOS_COV_DIR / "{sample}.intron_loci.bed",
        bam = BAM_DIR / "{sample}.sorted.bam",
    output:
        depth = NON_INTRONER_CIRCOS_COV_DIR / "{sample}.intron_depth.tsv",
    params:
        mapq = MAPQ_THRESHOLD,
    conda: "../envs/non_introner_introns.yaml"
    shell:
        """
        awk '!/^#/ && !/^contig/' {input.bed} | \
        samtools depth -Q {params.mapq} -a -J -b - {input.bam} > {output.depth}
        """


rule call_intron_presence_with_mating:
    """
    Call intron presence/absence for Circos-only with-mating intron loci.
    """
    input:
        bed = NON_INTRONER_CIRCOS_COV_DIR / "{sample}.intron_loci.bed",
        depth = NON_INTRONER_CIRCOS_COV_DIR / "{sample}.intron_depth.tsv",
    output:
        calls = NON_INTRONER_CIRCOS_COV_DIR / "{sample}.intron_coverage_calls.tsv",
    conda: "../envs/non_introner_introns.yaml"
    shell:
        """
        introner-caller --format simple {input.bed} {input.depth} {output.calls}
        """


rule make_intron_coverage_beds:
    """
    Create per-sample BED files of intron loci with 100bp flanking
    for the coverage caller. Only includes genes from reference catalog.
    """
    input:
        gtf = lambda wildcards: get_gtf(wildcards.sample),
        catalog = NON_INTRONER_DIR / "reference_intron_catalog.tsv",
        assembly = ASSEMBLIES_DIR / "{sample}.vg_paths.fa",
    output:
        bed = NON_INTRONER_COV_DIR / "{sample}.intron_loci.bed",
    params:
        sample = "{sample}",
    conda: "../envs/non_introner_introns.yaml"
    shell:
        """
        mkdir -p {NON_INTRONER_COV_DIR} && \
        python {NON_INTRONER_SCRIPTS}/make_intron_coverage_beds.py \
            --gtf {input.gtf} \
            --reference_catalog {input.catalog} \
            --assembly {input.assembly} \
            --sample {params.sample} \
            --output {output.bed}
        """


rule calculate_intron_locus_depth:
    """
    Calculate read depth at each intron locus using samtools depth.
    Reuses existing BAMs from the introner coverage calling pipeline.
    """
    input:
        bed = NON_INTRONER_COV_DIR / "{sample}.intron_loci.bed",
        bam = BAM_DIR / "{sample}.sorted.bam",
    output:
        depth = NON_INTRONER_COV_DIR / "{sample}.intron_depth.tsv",
    params:
        mapq = MAPQ_THRESHOLD,
    conda: "../envs/non_introner_introns.yaml"
    shell:
        """
        awk '!/^#/ && !/^contig/' {input.bed} | \
        samtools depth -Q {params.mapq} -a -J -b - {input.bam} > {output.depth}
        """


rule call_intron_presence:
    """
    Call intron presence/absence based on coverage patterns using
    introner-caller with --format simple (6-column BED input).
    """
    input:
        bed = NON_INTRONER_COV_DIR / "{sample}.intron_loci.bed",
        depth = NON_INTRONER_COV_DIR / "{sample}.intron_depth.tsv",
    output:
        calls = NON_INTRONER_COV_DIR / "{sample}.intron_coverage_calls.tsv",
    conda: "../envs/non_introner_introns.yaml"
    shell:
        """
        introner-caller --format simple {input.bed} {input.depth} {output.calls}
        """


rule non_introner_intron_afs_with_mating_for_circos:
    """
    Build a Circos-only non-introner intron genotype matrix that retains
    CCMP1545 scaffold 2 mating-type region introns.

    The standard non_introner_intron_afs rule remains filtered for evolutionary
    analyses. This alternate matrix is only for complete chromosome density
    tracks in the Circos plot.
    """
    input:
        catalog = NON_INTRONER_CIRCOS_DIR / "reference_intron_catalog.with_mating.tsv",
        gtfs = [str(get_gtf(s)) for s in GROUP1_SAMPLES + GROUP2_SAMPLES],
        assemblies = expand(
            str(ASSEMBLIES_DIR / "{sample}.vg_paths.fa"),
            sample=GROUP1_SAMPLES + GROUP2_SAMPLES,
        ),
        introner_loci = GENOTYPING_DIR / "introner_loci.bed",
        coverage_calls = expand(
            str(NON_INTRONER_CIRCOS_COV_DIR / "{sample}.intron_coverage_calls.tsv"),
            sample=NON_REF_SAMPLES,
        ),
    output:
        matrix = NON_INTRONER_CIRCOS_DIR / "intron_genotype_matrix.with_mating.tsv",
        afs_plot = FIGURES_DIR / "genome_alignment" / "non_introner_intron_afs.with_mating.pdf",
        afs_png = FIGURES_DIR / "genome_alignment" / "non_introner_intron_afs.with_mating.png",
        summary = NON_INTRONER_CIRCOS_DIR / "summary_statistics.with_mating.txt",
    params:
        gtf_dir = str(ANNOTATIONS_DIR),
        assemblies_dir = str(ASSEMBLIES_DIR),
        group1 = ",".join(GROUP1_SAMPLES),
        group2 = ",".join(GROUP2_SAMPLES),
        reference = REFERENCE,
        coverage_dir = str(NON_INTRONER_CIRCOS_COV_DIR),
        # Override default {sample}.gtf lookup for samples with augmented GTFs
        gtf_overrides = ",".join(
            f"{s}={get_gtf(s)}" for s in AUGMENTED_GTF_SAMPLES
        ),
    conda: "../envs/non_introner_introns.yaml"
    shell:
        """
        mkdir -p {NON_INTRONER_CIRCOS_DIR} {FIGURES_DIR}/genome_alignment && \
        python {NON_INTRONER_SCRIPTS}/non_introner_intron_afs.py \
            --gtf_dir {params.gtf_dir} \
            --assemblies_dir {params.assemblies_dir} \
            --introner_loci {input.introner_loci} \
            --group1 {params.group1} \
            --group2 {params.group2} \
            --reference {params.reference} \
            --output_matrix {output.matrix} \
            --output_afs {output.afs_plot} \
            --output_summary {output.summary} \
            --coverage_dir {params.coverage_dir} \
            --gtf_overrides {params.gtf_overrides} \
            --include_mating_region \
            --keep_reference_singletons
        """


rule non_introner_intron_afs:
    """
    Build non-introner intron genotype matrix and allele frequency spectra.
    Uses minimap2 alignment for orthologous intron matching, validated by
    short-read coverage to eliminate assembly-derived artifacts.
    """
    input:
        gtfs = [str(get_gtf(s)) for s in GROUP1_SAMPLES + GROUP2_SAMPLES],
        assemblies = expand(
            str(ASSEMBLIES_DIR / "{sample}.vg_paths.fa"),
            sample=GROUP1_SAMPLES + GROUP2_SAMPLES,
        ),
        introner_loci = GENOTYPING_DIR / "introner_loci.bed",
        coverage_calls = expand(
            str(NON_INTRONER_COV_DIR / "{sample}.intron_coverage_calls.tsv"),
            sample=NON_REF_SAMPLES,
        ),
    output:
        matrix = NON_INTRONER_DIR / "intron_genotype_matrix.tsv",
        afs_plot = FIGURES_DIR / "non_introner_intron_afs.pdf",
        afs_png = FIGURES_DIR / "non_introner_intron_afs.png",
        summary = NON_INTRONER_DIR / "summary_statistics.txt",
    params:
        gtf_dir = str(ANNOTATIONS_DIR),
        assemblies_dir = str(ASSEMBLIES_DIR),
        group1 = ",".join(GROUP1_SAMPLES),
        group2 = ",".join(GROUP2_SAMPLES),
        reference = REFERENCE,
        coverage_dir = str(NON_INTRONER_COV_DIR),
        # Override default {sample}.gtf lookup for samples with augmented GTFs
        gtf_overrides = ",".join(
            f"{s}={get_gtf(s)}" for s in AUGMENTED_GTF_SAMPLES
        ),
    conda: "../envs/non_introner_introns.yaml"
    shell:
        """
        mkdir -p {NON_INTRONER_DIR} {FIGURES_DIR} && \
        python {NON_INTRONER_SCRIPTS}/non_introner_intron_afs.py \
            --gtf_dir {params.gtf_dir} \
            --assemblies_dir {params.assemblies_dir} \
            --introner_loci {input.introner_loci} \
            --group1 {params.group1} \
            --group2 {params.group2} \
            --reference {params.reference} \
            --output_matrix {output.matrix} \
            --output_afs {output.afs_plot} \
            --output_summary {output.summary} \
            --coverage_dir {params.coverage_dir} \
            --gtf_overrides {params.gtf_overrides}
        """
