# Between-gene R2C2 association reported in Results and Figure 5B.
rule fit_r2c2_between_gene_expression:
    input:
        script = EXPRESSION_SCRIPT_DIR / "r2c2_expression/fit_between_gene_expression.py",
        prepare = EXPRESSION_SCRIPT_DIR / "r2c2_expression/prepare_r2c2_expression.py",
        helper = EXPRESSION_SCRIPT_DIR / "negative_binomial_expression_pilot/fit_expression_pilot_models.py",
        data = R2C2_DATA / "r2c2/expression_model_data.tsv"
    output:
        coefficients = FUNCTIONAL_DIR / "r2c2_between_gene/model_coefficients.tsv",
        support = FUNCTIONAL_DIR / "r2c2_between_gene/model_support.tsv",
        summary = FUNCTIONAL_DIR / "r2c2_between_gene/summary.txt"
    params: outdir = FUNCTIONAL_DIR / "r2c2_between_gene"
    log: EXPRESSION_LOG_DIR / "fit_r2c2_between_gene_expression.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --input {input.data:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """
