# Exon-anchored canonical-intron orthology, regional consensus and body Dxy (1).

rule canonical_intron_body_dxy:
    input:
        script = PROJECT_ROOT / "scripts/evolution/canonical_intron_dxy/run_analysis.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        canonical = EVOLUTION_DIR / "non_introner_introns/intron_genotype_matrix.tsv",
        catalog = EVOLUTION_DIR / "non_introner_introns/reference_intron_catalog.tsv",
        gtfs = expand(ANNOTATIONS_DIR / "{sample}.gtf", sample=ALL_SAMPLES),
        assemblies = expand(ASSEMBLIES_DIR / "{sample}.vg_paths.fa", sample=ALL_SAMPLES),
        indices = expand(ASSEMBLIES_DIR / "{sample}.vg_paths.fa.fai", sample=ALL_SAMPLES),
        bams = expand(GENOTYPING_DIR / "coverage/bams/{sample}.sorted.bam", sample=[s for s in ALL_SAMPLES if s != "CCMP1545"]),
        bais = expand(GENOTYPING_DIR / "coverage/bams/{sample}.sorted.bam.bai", sample=[s for s in ALL_SAMPLES if s != "CCMP1545"]),
        qc = expand(EVOLUTION_DIR / "introner_body_consensus/{sample}.introner_body.qc.tsv", sample=ALL_SAMPLES),
        combined_qc = EVOLUTION_DIR / "introner_body_consensus/introner_body_consensus.qc.tsv",
        vcfs = expand(EVOLUTION_DIR / "introner_body_consensus/{sample}.introner_body.filtered.vcf.gz", sample=[s for s in ALL_SAMPLES if s != "CCMP1545"]),
        vcf_indices = expand(EVOLUTION_DIR / "introner_body_consensus/{sample}.introner_body.filtered.vcf.gz.tbi", sample=[s for s in ALL_SAMPLES if s != "CCMP1545"]),
        body = EVOLUTION_DIR / "diversity_metrics/shared_introner_body_dxy_200bp.tsv",
        flanks = EVOLUTION_DIR / "diversity_metrics/all_samples_diversity_metrics_200bp.tsv",
        regional_script = PROJECT_ROOT / "scripts/evolution/canonical_intron_dxy/check_introner_regional_consensus.py",
        qc_script = PROJECT_ROOT / "scripts/evolution/qc_introner_body_consensus.py"
    output:
        report = directory(EVOLUTION_DIR / "canonical_intron_dxy")
    params:
        outdir = EVOLUTION_DIR / "canonical_intron_dxy"
    threads: 4
    log:
        RESULTS / "logs/promoted_analyses/canonical_intron_body_dxy.log"
    conda: "../envs/analysis_tools.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --output {params.outdir:q} --stage all --workers {threads} > {log:q} 2>&1
        """
