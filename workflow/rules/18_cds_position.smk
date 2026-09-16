# Manuscript-relevant analysis promoted from analysis/.

rule cds_three_prime_position_test:
    input:
        script = PROJECT_ROOT / "scripts/evolution/frequency_turnover/analyze_three_prime_cds_position.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        canonical = EVOLUTION_DIR / "non_introner_introns/intron_genotype_matrix.tsv",
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        colors = PROJECT_ROOT / "master_figure_color_guide.tsv",
        color_script = PROJECT_ROOT / "scripts/figure_color_guide.py"
    output:
        report = directory(EVOLUTION_DIR / "frequency_turnover" / "cds_three_prime_position_test")
    params:
        outdir = EVOLUTION_DIR / "frequency_turnover" / "cds_three_prime_position_test"
    threads: 1
    log:
        RESULTS / "logs/promoted_analyses/cds_three_prime_position_test.log"
    conda: "../envs/analysis_tools.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --introner-matrix {input.matrix:q} --non-introner-matrix {input.canonical:q} --gtf {input.gtf:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

