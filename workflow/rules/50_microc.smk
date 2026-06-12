# ============================================================
# 50_microc.smk - Micro-C Processing
# ============================================================
#
# Project-native Micro-C read processing and contact matrix generation.
# The default config points to the smaller pilot dataset; switching to the
# larger CC15353 dataset only requires changing microc.active_dataset.
#
# Outputs:
# - results/microc/bams_raw/{sample}.aligned.bam
# - results/microc/bams_final/{sample}.valid.bam
# - results/microc/pairs/{sample}.valid.pairs.gz
# - results/microc/matrix/{sample}.mcool
# - results/microc/qc/{sample}.dedup_stats.txt
#
# ============================================================

import glob
from pathlib import Path


MICROC_CONFIG = config.get("microc", {})
MICROC_SAMPLES = MICROC_CONFIG.get("samples", [])
MICROC_ACTIVE_DATASET = config.get(
    "microc_dataset", MICROC_CONFIG.get("active_dataset", "pilot")
)
MICROC_DATASETS = MICROC_CONFIG.get("datasets", {})
MICROC_DATASET_CONFIG = MICROC_DATASETS.get(MICROC_ACTIVE_DATASET, {})

if not MICROC_SAMPLES:
    raise ValueError("config.yaml microc.samples must contain at least one sample")
if not MICROC_DATASET_CONFIG:
    raise ValueError(
        f"config.yaml microc.datasets does not define active dataset "
        f"'{MICROC_ACTIVE_DATASET}'"
    )

MICROC_FASTQ_DIR = Path(MICROC_DATASET_CONFIG["fastq_dir"])
MICROC_READ_PATTERN = MICROC_DATASET_CONFIG.get(
    "read_pattern", "{fastq_id}_S*_{read}_001.fastq.gz"
)
MICROC_FASTQ_IDS = MICROC_DATASET_CONFIG.get("sample_fastq_ids", {})

MICROC_ALIGN_CONFIG = MICROC_CONFIG.get("alignment", {})
MICROC_MATRIX_CONFIG = MICROC_CONFIG.get("matrix", {})
MICROC_QC_CONFIG = MICROC_CONFIG.get("qc", {})
MICROC_TRACK_CONFIG = MICROC_CONFIG.get("tracks", {})
MICROC_INTRONER_PROFILE_CONFIG = MICROC_CONFIG.get("introner_profiles", {})

MICROC_RAW_BAM_DIR = MICROC_DIR / "bams_raw"
MICROC_FINAL_BAM_DIR = MICROC_DIR / "bams_final"
MICROC_PAIRS_DIR = MICROC_DIR / "pairs"
MICROC_TEMP_DIR = MICROC_DIR / "temp"
MICROC_QC_DIR = MICROC_DIR / "qc"
MICROC_MATRIX_DIR = MICROC_DIR / "matrix"
MICROC_TRACK_DIR = MICROC_DIR / "tracks"
MICROC_REF_DIR = MICROC_DIR / "references"
MICROC_LOG_DIR = MICROC_DIR / "logs"
MICROC_HEATMAP_DIR = MICROC_QC_DIR / "heatmaps"
MICROC_DECAY_DIR = MICROC_QC_DIR / "distance_decay"
MICROC_INSULATION_DIR = MICROC_QC_DIR / "insulation"
MICROC_INTRONER_PROFILE_DIR = MICROC_QC_DIR / "introner_profiles"
MICROC_CHROM_TRACK_DIR = MICROC_QC_DIR / "chromosome_tracks"

MICROC_HEATMAP_RES = MICROC_QC_CONFIG.get("heatmap_resolution", 5000)
MICROC_DETAIL_RES = MICROC_QC_CONFIG.get("detail_resolution", 1000)
MICROC_DECAY_RES = MICROC_QC_CONFIG.get("decay_resolution", 1000)
MICROC_INSULATION_RES = MICROC_QC_CONFIG.get("insulation_resolution", 1000)
MICROC_INSULATION_WINDOW = MICROC_QC_CONFIG.get("insulation_window", 25000)
MICROC_INSULATION_THRESHOLD_METHOD = MICROC_QC_CONFIG.get(
    "insulation_boundary_threshold_method", "li"
)
MICROC_INTRONER_BOUNDARY_FLANK = MICROC_QC_CONFIG.get("introner_boundary_flank", 5000)
MICROC_INTRONER_BOUNDARY_PERMUTATIONS = MICROC_QC_CONFIG.get(
    "introner_boundary_permutations", 1000
)
MICROC_INTRONER_INSULATION_PERMUTATIONS = MICROC_QC_CONFIG.get(
    "introner_insulation_permutations", 10000
)
MICROC_INTRONER_BOUNDARY_SEED = MICROC_QC_CONFIG.get(
    "introner_boundary_random_seed", 1545
)
MICROC_INTRONER_DOMAIN_POSITION_MIN_CHROM_SIZE = MICROC_QC_CONFIG.get(
    "introner_domain_position_min_chrom_size", 100000
)
MICROC_INTRONER_BOUNDARY_EXCLUDED_CHROMS = MICROC_QC_CONFIG.get(
    "introner_boundary_excluded_chroms", {}
)
MICROC_CHROM_TRACK_SAMPLES = MICROC_QC_CONFIG.get(
    "chromosome_track_samples", [REFERENCE]
)
MICROC_CHROM_TRACK_CHROMS = [
    str(chrom) for chrom in MICROC_QC_CONFIG.get("chromosome_track_chroms", [1, 3, 4, 5])
]
MICROC_CHROM_TRACK_MAX_DISTANCE = MICROC_QC_CONFIG.get(
    "chromosome_track_max_distance", 250000
)
MICROC_CHROM_TRACK_INTRONER_BIN = MICROC_QC_CONFIG.get(
    "chromosome_track_introner_bin_bp", 25000
)
MICROC_TRACK_MAPQ = MICROC_TRACK_CONFIG.get("mapq", 40)
MICROC_TRACK_EXCLUDE_FLAGS = MICROC_TRACK_CONFIG.get("exclude_flags", 3844)
MICROC_TRACK_THREADS = MICROC_TRACK_CONFIG.get("threads", 8)
MICROC_INTRONER_PROFILE_SAMPLES = MICROC_INTRONER_PROFILE_CONFIG.get(
    "samples", [REFERENCE]
)
MICROC_INTRONER_PROFILE_WINDOW = MICROC_INTRONER_PROFILE_CONFIG.get("window_bp", 1000)
MICROC_INTRONER_PROFILE_BIN = MICROC_INTRONER_PROFILE_CONFIG.get("bin_bp", 10)
MICROC_INTRONER_PROFILE_ABSENT_MAX = MICROC_INTRONER_PROFILE_CONFIG.get(
    "absent_max_length", 200
)
MICROC_NORMAL_INTRON_EXCLUSION_WINDOW = MICROC_INTRONER_PROFILE_CONFIG.get(
    "normal_intron_exclusion_window_bp", 1000
)


def get_microc_normal_intron_file(wildcards):
    """Return normal-intron intervals in the coordinate system of the sample."""
    if wildcards.sample == REFERENCE:
        return EVOLUTION_DIR / "non_introner_introns" / "reference_intron_catalog.tsv"
    standard_path = (
        EVOLUTION_DIR
        / "non_introner_introns"
        / "coverage"
        / f"{wildcards.sample}.intron_loci.bed"
    )
    circos_path = (
        EVOLUTION_DIR
        / "non_introner_introns"
        / "circos"
        / "coverage"
        / f"{wildcards.sample}.intron_loci.bed"
    )
    return standard_path if standard_path.exists() else circos_path


def get_microc_introner_boundary_excluded_chroms(wildcards):
    excluded = MICROC_INTRONER_BOUNDARY_EXCLUDED_CHROMS.get(wildcards.sample, [])
    if isinstance(excluded, str):
        excluded = [excluded]
    return ",".join(excluded)


def get_microc_fastq(wildcards, read):
    """Return the unique FASTQ for a Micro-C sample/read in the active dataset."""
    if wildcards.sample not in MICROC_FASTQ_IDS:
        raise ValueError(
            f"No FASTQ ID configured for Micro-C sample '{wildcards.sample}' "
            f"in dataset '{MICROC_ACTIVE_DATASET}'"
        )

    fastq_id = MICROC_FASTQ_IDS[wildcards.sample]
    pattern = MICROC_FASTQ_DIR / MICROC_READ_PATTERN.format(
        sample=wildcards.sample,
        fastq_id=fastq_id,
        read=read,
    )
    found = sorted(glob.glob(str(pattern)))

    if not found:
        raise ValueError(
            f"Could not find Micro-C {read} FASTQ for {wildcards.sample}; "
            f"searched pattern: {pattern}"
        )
    if len(found) > 1:
        raise ValueError(
            f"Found multiple Micro-C {read} FASTQs for {wildcards.sample}; "
            f"searched pattern: {pattern}; found: {found}"
        )

    return found[0]


rule microc_chrom_sizes:
    """
    Create chromosome sizes for Micro-C pair parsing and cooler binning.
    """
    input:
        fai = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.fai"
    output:
        chromsizes = MICROC_REF_DIR / "{sample}.chrom.sizes"
    shell:
        """
        mkdir -p {MICROC_REF_DIR}
        cut -f1,2 {input.fai} > {output.chromsizes}
        """


rule microc_link_reference:
    """
    Link the sample assembly into the Micro-C reference directory.

    bwa-mem2 writes index sidecars next to the reference. Keeping a Micro-C
    reference link avoids clobbering the classic BWA indexes in results/assemblies.
    """
    input:
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        fa = MICROC_REF_DIR / "{sample}.vg_paths.fa"
    shell:
        """
        mkdir -p {MICROC_REF_DIR}
        ln -sf {input.fa} {output.fa}
        """


rule microc_bwa_mem2_index:
    """
    Create a bwa-mem2 index for the sample-specific assembly.
    """
    input:
        fa = MICROC_REF_DIR / "{sample}.vg_paths.fa"
    output:
        idx0123 = MICROC_REF_DIR / "{sample}.vg_paths.fa.0123",
        amb = MICROC_REF_DIR / "{sample}.vg_paths.fa.amb",
        ann = MICROC_REF_DIR / "{sample}.vg_paths.fa.ann",
        bwt = MICROC_REF_DIR / "{sample}.vg_paths.fa.bwt.2bit.64",
        pac = MICROC_REF_DIR / "{sample}.vg_paths.fa.pac"
    log:
        MICROC_LOG_DIR / "bwa_mem2_index" / "{sample}.log"
    conda:
        "../envs/microc.yaml"
    shell:
        """
        mkdir -p {MICROC_LOG_DIR}/bwa_mem2_index
        bwa-mem2 index {input.fa} > {log} 2>&1
        """


rule microc_bwa_mem_raw:
    """
    Align Micro-C reads
    """
    input:
        ref = MICROC_REF_DIR / "{sample}.vg_paths.fa",
        idx = rules.microc_bwa_mem2_index.output,
        r1 = lambda wildcards: get_microc_fastq(wildcards, "R1"),
        r2 = lambda wildcards: get_microc_fastq(wildcards, "R2")
    output:
        bam = MICROC_RAW_BAM_DIR / "{sample}.aligned.bam"
    threads:
        MICROC_ALIGN_CONFIG.get("threads", 8)
        + MICROC_ALIGN_CONFIG.get("samtools_view_threads", 4)
    params:
        bwa_threads = MICROC_ALIGN_CONFIG.get("threads", 8),
        bwa_opts = MICROC_ALIGN_CONFIG.get("bwa_opts", "-5SP"),
        view_threads = MICROC_ALIGN_CONFIG.get("samtools_view_threads", 4)
    log:
        MICROC_LOG_DIR / "bwa_mem2" / "{sample}.log"
    conda:
        "../envs/microc.yaml"
    shell:
        """
        mkdir -p {MICROC_RAW_BAM_DIR} {MICROC_LOG_DIR}/bwa_mem2
        rm -f {output.bam}.tmp
        (
            bwa-mem2 mem {params.bwa_opts} -t {params.bwa_threads} {input.ref} {input.r1} {input.r2} \
            | samtools view -@ {params.view_threads} -bS - > {output.bam}.tmp

            samtools quickcheck {output.bam}.tmp
            mv {output.bam}.tmp {output.bam}
        ) 2> {log}
        """


rule microc_parse_bam:
    """
    Parse aligned BAM into pairsam format with Micro-C settings.
    """
    input:
        bam = MICROC_RAW_BAM_DIR / "{sample}.aligned.bam",
        chromsizes = MICROC_REF_DIR / "{sample}.chrom.sizes"
    output:
        pairsam = temp(MICROC_TEMP_DIR / "{sample}.parsed.pairsam.gz")
    threads:
        MICROC_ALIGN_CONFIG.get("threads", 8)
    params:
        mapq = MICROC_ALIGN_CONFIG.get("mapq", 40),
        walks_policy = MICROC_ALIGN_CONFIG.get("walks_policy", "5unique"),
        max_inter_align_gap = MICROC_ALIGN_CONFIG.get("max_inter_align_gap", 30)
    log:
        MICROC_LOG_DIR / "pairtools_parse" / "{sample}.log"
    conda:
        "../envs/microc.yaml"
    shell:
        """
        mkdir -p {MICROC_TEMP_DIR} {MICROC_LOG_DIR}/pairtools_parse
        pairtools parse \
            --min-mapq {params.mapq} \
            --walks-policy {params.walks_policy} \
            --max-inter-align-gap {params.max_inter_align_gap} \
            --chroms-path {input.chromsizes} \
            --output {output.pairsam} \
            --assembly {wildcards.sample} \
            --no-flip \
            {input.bam} \
            > {log} 2>&1
        """


rule microc_sort_pairs:
    """
    Sort parsed pairsam records.
    """
    input:
        pairsam = MICROC_TEMP_DIR / "{sample}.parsed.pairsam.gz"
    output:
        pairsam = temp(MICROC_TEMP_DIR / "{sample}.sorted.pairsam.gz")
    threads:
        MICROC_ALIGN_CONFIG.get("threads", 8)
    log:
        MICROC_LOG_DIR / "pairtools_sort" / "{sample}.log"
    conda:
        "../envs/microc.yaml"
    shell:
        """
        mkdir -p {MICROC_LOG_DIR}/pairtools_sort
        pairtools sort \
            --nproc {threads} \
            --output {output.pairsam} \
            {input.pairsam} \
            > {log} 2>&1
        """


rule microc_dedup_pairs:
    """
    Deduplicate sorted Micro-C pairs.
    """
    input:
        pairsam = MICROC_TEMP_DIR / "{sample}.sorted.pairsam.gz"
    output:
        pairsam = temp(MICROC_TEMP_DIR / "{sample}.dedup.pairsam.gz"),
        stats = MICROC_QC_DIR / "{sample}.dedup_stats.txt"
    threads:
        MICROC_ALIGN_CONFIG.get("threads", 8)
    log:
        MICROC_LOG_DIR / "pairtools_dedup" / "{sample}.log"
    conda:
        "../envs/microc.yaml"
    shell:
        """
        mkdir -p {MICROC_QC_DIR} {MICROC_LOG_DIR}/pairtools_dedup
        pairtools dedup \
            --mark-dups \
            --output-stats {output.stats} \
            --output {output.pairsam} \
            {input.pairsam} \
            > {log} 2>&1
        """


rule microc_split_bam_pairs:
    """
    Split deduplicated pairsam into final pairs and coordinate-sorted BAM.
    """
    input:
        pairsam = MICROC_TEMP_DIR / "{sample}.dedup.pairsam.gz"
    output:
        pairs = MICROC_PAIRS_DIR / "{sample}.valid.pairs.gz",
        bam = MICROC_FINAL_BAM_DIR / "{sample}.valid.bam",
        bai = MICROC_FINAL_BAM_DIR / "{sample}.valid.bam.bai"
    threads:
        MICROC_ALIGN_CONFIG.get("threads", 8)
    log:
        MICROC_LOG_DIR / "pairtools_split" / "{sample}.log"
    conda:
        "../envs/microc.yaml"
    shell:
        """
        mkdir -p {MICROC_PAIRS_DIR} {MICROC_FINAL_BAM_DIR} {MICROC_LOG_DIR}/pairtools_split
        (
            pairtools split \
                --output-pairs {output.pairs} \
                --output-sam - \
                {input.pairsam} \
            | samtools sort -@ {threads} -o {output.bam} -

            samtools index {output.bam}
        ) > {log} 2>&1
        """


rule microc_cload_pairs:
    """
    Create base-resolution cooler from valid Micro-C pairs.
    """
    input:
        pairs = MICROC_PAIRS_DIR / "{sample}.valid.pairs.gz",
        chromsizes = MICROC_REF_DIR / "{sample}.chrom.sizes"
    output:
        cool = temp(MICROC_MATRIX_DIR / "{sample}.base.cool")
    params:
        bin_size = MICROC_MATRIX_CONFIG.get("base_resolution", 200)
    log:
        MICROC_LOG_DIR / "cooler_cload" / "{sample}.log"
    conda:
        "../envs/microc.yaml"
    shell:
        """
        mkdir -p {MICROC_MATRIX_DIR} {MICROC_LOG_DIR}/cooler_cload
        cooler cload pairs \
            --chrom1 2 --pos1 3 \
            --chrom2 4 --pos2 5 \
            {input.chromsizes}:{params.bin_size} \
            {input.pairs} \
            {output.cool} \
            > {log} 2>&1
        """


rule microc_zoomify:
    """
    Create balanced multi-resolution contact matrix.
    """
    input:
        cool = MICROC_MATRIX_DIR / "{sample}.base.cool"
    output:
        mcool = MICROC_MATRIX_DIR / "{sample}.mcool"
    threads:
        MICROC_MATRIX_CONFIG.get("threads", 8)
    params:
        resolutions = ",".join(map(str, MICROC_MATRIX_CONFIG.get(
            "resolutions", [200, 400, 1000, 2000, 5000, 10000, 25000, 50000]
        )))
    log:
        MICROC_LOG_DIR / "cooler_zoomify" / "{sample}.log"
    conda:
        "../envs/microc.yaml"
    shell:
        """
        mkdir -p {MICROC_LOG_DIR}/cooler_zoomify
        cooler zoomify \
            --balance \
            --nproc {threads} \
            --resolutions {params.resolutions} \
            --out {output.mcool} \
            {input.cool} \
            > {log} 2>&1
        """


rule microc_contact_heatmaps:
    """
    Create global and per-contig Micro-C contact heatmaps for visual QC.
    """
    input:
        mcool = MICROC_MATRIX_DIR / "{sample}.mcool",
        script = PROJECT_ROOT / "scripts" / "microc" / "plot_contact_heatmaps.py"
    output:
        global_png = MICROC_HEATMAP_DIR / f"{{sample}}.global_{MICROC_HEATMAP_RES}bp.png",
        global_pdf = MICROC_HEATMAP_DIR / f"{{sample}}.global_{MICROC_HEATMAP_RES}bp.pdf",
        manifest = MICROC_HEATMAP_DIR / "{sample}.heatmaps.tsv"
    params:
        outdir = MICROC_HEATMAP_DIR,
        global_resolution = MICROC_HEATMAP_RES,
        detail_resolution = MICROC_DETAIL_RES,
        max_detail_chroms = MICROC_QC_CONFIG.get("detail_max_chroms", 12),
        min_detail_size = MICROC_QC_CONFIG.get("detail_min_chrom_size", 100000)
    log:
        MICROC_LOG_DIR / "microc_contact_heatmaps" / "{sample}.log"
    conda:
        "../envs/microc_qc.yaml"
    shell:
        """
        mkdir -p {MICROC_HEATMAP_DIR} {MICROC_LOG_DIR}/microc_contact_heatmaps
        python {input.script} \
            --mcool {input.mcool} \
            --sample {wildcards.sample} \
            --global-resolution {params.global_resolution} \
            --detail-resolution {params.detail_resolution} \
            --max-detail-chroms {params.max_detail_chroms} \
            --min-detail-size {params.min_detail_size} \
            --outdir {params.outdir} \
            > {log} 2>&1
        """


rule microc_distance_decay:
    """
    Estimate and plot cis contact decay P(s) from Micro-C matrices.
    """
    input:
        mcool = MICROC_MATRIX_DIR / "{sample}.mcool",
        script = PROJECT_ROOT / "scripts" / "microc" / "plot_distance_decay.py"
    output:
        by_diagonal = MICROC_DECAY_DIR / f"{{sample}}.{MICROC_DECAY_RES}bp.by_diagonal.tsv",
        binned = MICROC_DECAY_DIR / f"{{sample}}.{MICROC_DECAY_RES}bp.binned.tsv",
        plot_png = MICROC_DECAY_DIR / f"{{sample}}.{MICROC_DECAY_RES}bp.png",
        plot_pdf = MICROC_DECAY_DIR / f"{{sample}}.{MICROC_DECAY_RES}bp.pdf"
    params:
        out_prefix = lambda wildcards: MICROC_DECAY_DIR / f"{wildcards.sample}.{MICROC_DECAY_RES}bp",
        resolution = MICROC_DECAY_RES,
        max_distance = MICROC_QC_CONFIG.get("decay_max_distance", 2000000),
        log_bins = MICROC_QC_CONFIG.get("decay_log_bins", 80)
    log:
        MICROC_LOG_DIR / "microc_distance_decay" / "{sample}.log"
    conda:
        "../envs/microc_qc.yaml"
    shell:
        """
        mkdir -p {MICROC_DECAY_DIR} {MICROC_LOG_DIR}/microc_distance_decay
        python {input.script} \
            --mcool {input.mcool} \
            --sample {wildcards.sample} \
            --resolution {params.resolution} \
            --max-distance {params.max_distance} \
            --log-bins {params.log_bins} \
            --out-prefix {params.out_prefix} \
            > {log} 2>&1
        """


rule microc_insulation_domains:
    """
    Calculate insulation tracks, boundary BEDs, and domain BEDs from Micro-C matrices.
    """
    input:
        mcool = MICROC_MATRIX_DIR / "{sample}.mcool",
        script = PROJECT_ROOT / "scripts" / "microc" / "call_insulation_domains.py"
    output:
        insulation = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.insulation.tsv",
        boundaries = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.boundaries.bed",
        candidate_boundaries = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.candidate_boundaries.bed",
        cooltools_boundaries = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.cooltools_boundaries.bed",
        li_boundaries = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.li_boundaries.bed",
        otsu_boundaries = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.otsu_boundaries.bed",
        boundary_thresholds = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.boundary_thresholds.tsv",
        domains = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.domains.bed",
        plot_png = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.insulation.png",
        plot_pdf = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.insulation.pdf"
    params:
        out_prefix = lambda wildcards: MICROC_INSULATION_DIR / f"{wildcards.sample}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp",
        resolution = MICROC_INSULATION_RES,
        window = MICROC_INSULATION_WINDOW,
        ignore_diags = MICROC_QC_CONFIG.get("insulation_ignore_diags", 2),
        min_boundary_strength = MICROC_QC_CONFIG.get("insulation_min_boundary_strength", 0.0),
        threshold_method = MICROC_INSULATION_THRESHOLD_METHOD,
        max_plot_chroms = MICROC_QC_CONFIG.get("insulation_max_plot_chroms", 12),
        min_plot_size = MICROC_QC_CONFIG.get("insulation_min_plot_chrom_size", 100000)
    log:
        MICROC_LOG_DIR / "microc_insulation_domains" / "{sample}.log"
    conda:
        "../envs/microc_qc.yaml"
    shell:
        """
        mkdir -p {MICROC_INSULATION_DIR} {MICROC_LOG_DIR}/microc_insulation_domains
        python {input.script} \
            --mcool {input.mcool} \
            --sample {wildcards.sample} \
            --resolution {params.resolution} \
            --window {params.window} \
            --ignore-diags {params.ignore_diags} \
            --min-boundary-strength {params.min_boundary_strength} \
            --boundary-threshold-method {params.threshold_method} \
            --max-plot-chroms {params.max_plot_chroms} \
            --min-plot-size {params.min_plot_size} \
            --out-prefix {params.out_prefix} \
            > {log} 2>&1
        """


rule microc_introner_insulation_enrichment:
    """
    Test present-introner midpoint insulation against gene-conditioned permutations.
    """
    input:
        genotype = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        chromsizes = MICROC_REF_DIR / "{sample}.chrom.sizes",
        insulation = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.insulation.tsv",
        genes = lambda wildcards: get_gtf(wildcards.sample),
        script = PROJECT_ROOT / "scripts" / "microc" / "analyze_introner_boundary_enrichment.py"
    output:
        summary = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.introner_insulation_enrichment.tsv",
        chrom_summary = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.introner_insulation_by_chrom.tsv",
        loci = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.introner_insulation_loci.tsv",
        null = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.introner_insulation_null.tsv",
        plot_png = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.introner_insulation_enrichment.png",
        plot_pdf = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.introner_insulation_enrichment.pdf"
    params:
        excluded_chroms = get_microc_introner_boundary_excluded_chroms,
        permutations = MICROC_INTRONER_INSULATION_PERMUTATIONS,
        seed = MICROC_INTRONER_BOUNDARY_SEED
    log:
        MICROC_LOG_DIR / "microc_introner_insulation_enrichment" / "{sample}.log"
    conda:
        "../envs/microc_qc.yaml"
    shell:
        """
        mkdir -p {MICROC_INSULATION_DIR} {MICROC_LOG_DIR}/microc_introner_insulation_enrichment
        python {input.script} \
            --genotype-matrix {input.genotype} \
            --sample {wildcards.sample} \
            --chrom-sizes {input.chromsizes} \
            --insulation {input.insulation} \
            --genes-gtf {input.genes} \
            --exclude-chroms "{params.excluded_chroms}" \
            --window {MICROC_INSULATION_WINDOW} \
            --permutations {params.permutations} \
            --seed {params.seed} \
            --summary-tsv {output.summary} \
            --chrom-summary-tsv {output.chrom_summary} \
            --loci-tsv {output.loci} \
            --null-tsv {output.null} \
            --png {output.plot_png} \
            --pdf {output.plot_pdf} \
            > {log} 2>&1
        """


rule microc_chromosome_tad_tracks:
    """
    Plot chromosome-scale Micro-C diagonal, insulation, boundaries, and introner bins.
    """
    input:
        mcool = MICROC_MATRIX_DIR / "{sample}.mcool",
        chromsizes = MICROC_REF_DIR / "{sample}.chrom.sizes",
        genotype = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        genes = lambda wildcards: get_gtf(wildcards.sample),
        insulation = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.insulation.tsv",
        boundaries = MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.boundaries.bed",
        script = PROJECT_ROOT / "scripts" / "microc" / "plot_chromosome_tad_tracks.py"
    output:
        manifest = MICROC_CHROM_TRACK_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.chromosome_tad_tracks.manifest.tsv"
    params:
        out_prefix = lambda wildcards: MICROC_CHROM_TRACK_DIR / f"{wildcards.sample}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.chromosome_tad_tracks",
        chroms = ",".join(MICROC_CHROM_TRACK_CHROMS),
        max_distance = MICROC_CHROM_TRACK_MAX_DISTANCE,
        introner_bin = MICROC_CHROM_TRACK_INTRONER_BIN,
        absent_max = MICROC_INTRONER_PROFILE_ABSENT_MAX
    log:
        MICROC_LOG_DIR / "microc_chromosome_tad_tracks" / "{sample}.log"
    conda:
        "../envs/microc_qc.yaml"
    shell:
        """
        mkdir -p {MICROC_CHROM_TRACK_DIR} {MICROC_LOG_DIR}/microc_chromosome_tad_tracks
        python {input.script} \
            --mcool {input.mcool} \
            --sample {wildcards.sample} \
            --chrom-sizes {input.chromsizes} \
            --genotype-matrix {input.genotype} \
            --genes-gtf {input.genes} \
            --insulation {input.insulation} \
            --boundaries {input.boundaries} \
            --chroms {params.chroms} \
            --resolution {MICROC_INSULATION_RES} \
            --window {MICROC_INSULATION_WINDOW} \
            --max-distance {params.max_distance} \
            --introner-bin-bp {params.introner_bin} \
            --absent-max-length {params.absent_max} \
            --out-prefix {params.out_prefix} \
            > {log} 2>&1
        """


rule microc_track_scaling:
    """
    Count usable alignments and calculate the CPM scale factor for Micro-C tracks.
    """
    input:
        bam = MICROC_FINAL_BAM_DIR / "{sample}.valid.bam"
    output:
        stats = MICROC_TRACK_DIR / "{sample}.track_scaling.tsv"
    threads:
        MICROC_TRACK_THREADS
    params:
        mapq = MICROC_TRACK_MAPQ,
        exclude_flags = MICROC_TRACK_EXCLUDE_FLAGS
    log:
        MICROC_LOG_DIR / "microc_tracks" / "{sample}.scaling.log"
    conda:
        "../envs/microc_tracks.yaml"
    shell:
        """
        mkdir -p {MICROC_TRACK_DIR} {MICROC_LOG_DIR}/microc_tracks
        count=$(samtools view \
            -@ {threads} \
            -c \
            -q {params.mapq} \
            -F {params.exclude_flags} \
            {input.bam})
        scale=$(python -c 'import sys; c=int(sys.argv[1]); print(1000000.0 / c if c else 0)' "$count")
        printf "sample\tusable_alignments\tscale_cpm\tmapq\texclude_flags\n" > {output.stats}
        printf "{wildcards.sample}\t%s\t%s\t{params.mapq}\t{params.exclude_flags}\n" "$count" "$scale" >> {output.stats}
        printf "usable_alignments=%s\nscale_cpm=%s\n" "$count" "$scale" > {log}
        """


rule microc_anchor_track:
    """
    Build CPM-normalized 5-prime anchor tracks from valid Micro-C alignments.
    """
    input:
        bam = MICROC_FINAL_BAM_DIR / "{sample}.valid.bam",
        chromsizes = MICROC_REF_DIR / "{sample}.chrom.sizes",
        stats = MICROC_TRACK_DIR / "{sample}.track_scaling.tsv"
    output:
        bedgraph = MICROC_TRACK_DIR / "{sample}.anchors.cpm.bedgraph.gz",
        bigwig = MICROC_TRACK_DIR / "{sample}.anchors.cpm.bw"
    threads:
        MICROC_TRACK_THREADS
    params:
        mapq = MICROC_TRACK_MAPQ,
        exclude_flags = MICROC_TRACK_EXCLUDE_FLAGS
    log:
        MICROC_LOG_DIR / "microc_tracks" / "{sample}.anchors.log"
    conda:
        "../envs/microc_tracks.yaml"
    shell:
        """
        mkdir -p {MICROC_TRACK_DIR} {MICROC_LOG_DIR}/microc_tracks
        tmp_bg={output.bedgraph}.tmp
        rm -f "$tmp_bg" {output.bigwig}
        scale=$(awk 'NR == 2 {{print $3}}' {input.stats})
        samtools view \
            -@ {threads} \
            -b \
            -q {params.mapq} \
            -F {params.exclude_flags} \
            {input.bam} \
        | bedtools genomecov \
            -ibam stdin \
            -5 \
            -bg \
            -scale "$scale" \
            > "$tmp_bg"
        bedGraphToBigWig "$tmp_bg" {input.chromsizes} {output.bigwig}
        gzip -c "$tmp_bg" > {output.bedgraph}
        rm -f "$tmp_bg"
        """


rule microc_coverage_track:
    """
    Build CPM-normalized read coverage tracks from valid Micro-C alignments.
    """
    input:
        bam = MICROC_FINAL_BAM_DIR / "{sample}.valid.bam",
        chromsizes = MICROC_REF_DIR / "{sample}.chrom.sizes",
        stats = MICROC_TRACK_DIR / "{sample}.track_scaling.tsv"
    output:
        bedgraph = MICROC_TRACK_DIR / "{sample}.coverage.cpm.bedgraph.gz",
        bigwig = MICROC_TRACK_DIR / "{sample}.coverage.cpm.bw"
    threads:
        MICROC_TRACK_THREADS
    params:
        mapq = MICROC_TRACK_MAPQ,
        exclude_flags = MICROC_TRACK_EXCLUDE_FLAGS
    log:
        MICROC_LOG_DIR / "microc_tracks" / "{sample}.coverage.log"
    conda:
        "../envs/microc_tracks.yaml"
    shell:
        """
        mkdir -p {MICROC_TRACK_DIR} {MICROC_LOG_DIR}/microc_tracks
        tmp_bg={output.bedgraph}.tmp
        rm -f "$tmp_bg" {output.bigwig}
        scale=$(awk 'NR == 2 {{print $3}}' {input.stats})
        samtools view \
            -@ {threads} \
            -b \
            -q {params.mapq} \
            -F {params.exclude_flags} \
            {input.bam} \
        | bedtools genomecov \
            -ibam stdin \
            -bg \
            -scale "$scale" \
            > "$tmp_bg"
        bedGraphToBigWig "$tmp_bg" {input.chromsizes} {output.bigwig}
        gzip -c "$tmp_bg" > {output.bedgraph}
        rm -f "$tmp_bg"
        """


rule microc_introner_anchor_profile:
    """
    Aggregate Micro-C anchor signal around genotype-matrix introner boundaries.
    """
    input:
        genotype = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        anchors = MICROC_TRACK_DIR / "{sample}.anchors.cpm.bedgraph.gz",
        chromsizes = MICROC_REF_DIR / "{sample}.chrom.sizes",
        normal_introns = get_microc_normal_intron_file,
        script = PROJECT_ROOT / "scripts" / "microc" / "plot_introner_anchor_profiles.py"
    output:
        profile = MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_anchor_boundary_profile.tsv",
        intron_comparison_profile = MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_vs_normal_intron_anchor_boundary_profile.tsv",
        summary = MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_anchor_boundary_summary.tsv",
        bed = MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_states.bed",
        bedgz = MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_states.bed.gz",
        tbi = MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_states.bed.gz.tbi",
        png = MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_anchor_boundary_profile.png",
        pdf = MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_anchor_boundary_profile.pdf",
        intron_comparison_png = MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_vs_normal_intron_anchor_boundary_profile.png",
        intron_comparison_pdf = MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_vs_normal_intron_anchor_boundary_profile.pdf"
    params:
        window = MICROC_INTRONER_PROFILE_WINDOW,
        bin = MICROC_INTRONER_PROFILE_BIN,
        absent_max = MICROC_INTRONER_PROFILE_ABSENT_MAX,
        intron_exclusion_window = MICROC_NORMAL_INTRON_EXCLUSION_WINDOW
    log:
        MICROC_LOG_DIR / "microc_introner_profiles" / "{sample}.log"
    conda:
        "../envs/microc_qc.yaml"
    shell:
        """
        mkdir -p {MICROC_INTRONER_PROFILE_DIR} {MICROC_LOG_DIR}/microc_introner_profiles
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        mkdir -p "$MPLCONFIGDIR"
        python {input.script} \
            --genotype-matrix {input.genotype} \
            --sample {wildcards.sample} \
            --anchor-bedgraph {input.anchors} \
            --chrom-sizes {input.chromsizes} \
            --normal-introns {input.normal_introns} \
            --window-bp {params.window} \
            --bin-bp {params.bin} \
            --absent-max-length {params.absent_max} \
            --normal-intron-exclusion-window-bp {params.intron_exclusion_window} \
            --profile-tsv {output.profile} \
            --intron-comparison-profile-tsv {output.intron_comparison_profile} \
            --summary-tsv {output.summary} \
            --bed {output.bed} \
            --png {output.png} \
            --pdf {output.pdf} \
            --intron-comparison-png {output.intron_comparison_png} \
            --intron-comparison-pdf {output.intron_comparison_pdf} \
            > {log} 2>&1
        bgzip -f -c {output.bed} > {output.bedgz}
        tabix -f -p bed {output.bedgz}
        """


rule microc_processing:
    """
    Target: process Micro-C reads through deduplicated pairs and BAMs.
    """
    input:
        expand(MICROC_PAIRS_DIR / "{sample}.valid.pairs.gz", sample=MICROC_SAMPLES),
        expand(MICROC_FINAL_BAM_DIR / "{sample}.valid.bam", sample=MICROC_SAMPLES),
        expand(MICROC_FINAL_BAM_DIR / "{sample}.valid.bam.bai", sample=MICROC_SAMPLES),
        expand(MICROC_QC_DIR / "{sample}.dedup_stats.txt", sample=MICROC_SAMPLES)


rule microc_matrices:
    """
    Target: build Micro-C multi-resolution contact matrices.
    """
    input:
        expand(MICROC_MATRIX_DIR / "{sample}.mcool", sample=MICROC_SAMPLES)


rule microc_visual_qc:
    """
    Target: create global and per-contig Micro-C heatmaps.
    """
    input:
        expand(MICROC_HEATMAP_DIR / f"{{sample}}.global_{MICROC_HEATMAP_RES}bp.png", sample=MICROC_SAMPLES),
        expand(MICROC_HEATMAP_DIR / f"{{sample}}.global_{MICROC_HEATMAP_RES}bp.pdf", sample=MICROC_SAMPLES),
        expand(MICROC_HEATMAP_DIR / "{sample}.heatmaps.tsv", sample=MICROC_SAMPLES)


rule microc_decay_qc:
    """
    Target: calculate distance-decay/P(s) QC for Micro-C matrices.
    """
    input:
        expand(MICROC_DECAY_DIR / f"{{sample}}.{MICROC_DECAY_RES}bp.by_diagonal.tsv", sample=MICROC_SAMPLES),
        expand(MICROC_DECAY_DIR / f"{{sample}}.{MICROC_DECAY_RES}bp.binned.tsv", sample=MICROC_SAMPLES),
        expand(MICROC_DECAY_DIR / f"{{sample}}.{MICROC_DECAY_RES}bp.png", sample=MICROC_SAMPLES),
        expand(MICROC_DECAY_DIR / f"{{sample}}.{MICROC_DECAY_RES}bp.pdf", sample=MICROC_SAMPLES)


rule microc_insulation_qc:
    """
    Target: calculate insulation, boundary, and domain QC for Micro-C matrices.
    """
    input:
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.insulation.tsv", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.boundaries.bed", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.candidate_boundaries.bed", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.cooltools_boundaries.bed", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.li_boundaries.bed", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.otsu_boundaries.bed", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.boundary_thresholds.tsv", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.domains.bed", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.insulation.png", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.insulation.pdf", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.introner_insulation_enrichment.tsv", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.introner_insulation_by_chrom.tsv", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.introner_insulation_loci.tsv", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.introner_insulation_null.tsv", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.introner_insulation_enrichment.png", sample=MICROC_SAMPLES),
        expand(MICROC_INSULATION_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.introner_insulation_enrichment.pdf", sample=MICROC_SAMPLES)


rule microc_chromosome_tad_track_figures:
    """
    Target: plot selected chromosome tracks for Micro-C/TAD/introner interpretation.
    """
    input:
        expand(
            MICROC_CHROM_TRACK_DIR / f"{{sample}}.{MICROC_INSULATION_RES}bp_{MICROC_INSULATION_WINDOW}bp.chromosome_tad_tracks.manifest.tsv",
            sample=MICROC_CHROM_TRACK_SAMPLES,
        )


rule microc_qc:
    """
    Target: complete Micro-C matrix QC.
    """
    input:
        rules.microc_visual_qc.input,
        rules.microc_decay_qc.input,
        rules.microc_insulation_qc.input,
        rules.microc_chromosome_tad_track_figures.input


rule microc_nucleosome_tracks:
    """
    Target: build browser-ready Micro-C nucleosome-scale coverage and anchor tracks.
    """
    input:
        expand(MICROC_TRACK_DIR / "{sample}.track_scaling.tsv", sample=MICROC_SAMPLES),
        expand(MICROC_TRACK_DIR / "{sample}.anchors.cpm.bedgraph.gz", sample=MICROC_SAMPLES),
        expand(MICROC_TRACK_DIR / "{sample}.anchors.cpm.bw", sample=MICROC_SAMPLES),
        expand(MICROC_TRACK_DIR / "{sample}.coverage.cpm.bedgraph.gz", sample=MICROC_SAMPLES),
        expand(MICROC_TRACK_DIR / "{sample}.coverage.cpm.bw", sample=MICROC_SAMPLES)


rule microc_introner_profiles:
    """
    Target: plot aggregate Micro-C anchor signal around introner boundary classes.
    """
    input:
        expand(
            MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_anchor_boundary_profile.tsv",
            sample=MICROC_INTRONER_PROFILE_SAMPLES,
        ),
        expand(
            MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_vs_normal_intron_anchor_boundary_profile.tsv",
            sample=MICROC_INTRONER_PROFILE_SAMPLES,
        ),
        expand(
            MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_anchor_boundary_summary.tsv",
            sample=MICROC_INTRONER_PROFILE_SAMPLES,
        ),
        expand(
            MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_states.bed.gz",
            sample=MICROC_INTRONER_PROFILE_SAMPLES,
        ),
        expand(
            MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_states.bed.gz.tbi",
            sample=MICROC_INTRONER_PROFILE_SAMPLES,
        ),
        expand(
            MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_anchor_boundary_profile.png",
            sample=MICROC_INTRONER_PROFILE_SAMPLES,
        ),
        expand(
            MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_anchor_boundary_profile.pdf",
            sample=MICROC_INTRONER_PROFILE_SAMPLES,
        ),
        expand(
            MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_vs_normal_intron_anchor_boundary_profile.png",
            sample=MICROC_INTRONER_PROFILE_SAMPLES,
        ),
        expand(
            MICROC_INTRONER_PROFILE_DIR / "{sample}.introner_vs_normal_intron_anchor_boundary_profile.pdf",
            sample=MICROC_INTRONER_PROFILE_SAMPLES,
        )


rule microc_all:
    """
    Target: complete Micro-C processing, matrix generation, matrix QC, and signal tracks.
    """
    input:
        rules.microc_processing.input,
        rules.microc_matrices.input,
        rules.microc_qc.input,
        rules.microc_nucleosome_tracks.input,
        rules.microc_introner_profiles.input
