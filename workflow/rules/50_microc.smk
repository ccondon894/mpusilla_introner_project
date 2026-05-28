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

MICROC_RAW_BAM_DIR = MICROC_DIR / "bams_raw"
MICROC_FINAL_BAM_DIR = MICROC_DIR / "bams_final"
MICROC_PAIRS_DIR = MICROC_DIR / "pairs"
MICROC_TEMP_DIR = MICROC_DIR / "temp"
MICROC_QC_DIR = MICROC_DIR / "qc"
MICROC_MATRIX_DIR = MICROC_DIR / "matrix"
MICROC_REF_DIR = MICROC_DIR / "references"
MICROC_LOG_DIR = MICROC_DIR / "logs"


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


rule microc_all:
    """
    Target: complete Micro-C processing and matrix generation.
    """
    input:
        rules.microc_processing.input,
        rules.microc_matrices.input
