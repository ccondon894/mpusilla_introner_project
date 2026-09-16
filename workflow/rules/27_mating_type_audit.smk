# Alignment evidence supporting the operational RCC1749 MT interval (48).

rule rcc1749_mating_type_alignment_audit:
    input:
        script = PROJECT_ROOT / "scripts/genome_alignment/audit_mt_alignment.py",
        coords = GENOME_ALIGNMENT_DIR / "mummer/CCMP1545_vs_RCC1749.coords",
        gtf = ANNOTATIONS_DIR / "RCC1749.gtf",
        colors = PROJECT_ROOT / "master_figure_color_guide.tsv",
        helper = PROJECT_ROOT / "scripts/expression/negative_binomial_expression_pilot/prepare_expression_pilot_data.py"
    output:
        report = directory(GENOME_ALIGNMENT_DIR / "mating_type_audit")
    params:
        outdir = GENOME_ALIGNMENT_DIR / "mating_type_audit"
    threads: 1
    log:
        RESULTS / "logs/promoted_analyses/rcc1749_mating_type_alignment_audit.log"
    conda: "../envs/analysis_tools.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """
