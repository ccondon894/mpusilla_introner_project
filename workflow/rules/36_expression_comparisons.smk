# Between-gene R2C2 association reported in Results and Figure 5B.
rule fit_between_gene_expression:
    input:
        script = EXPRESSION_SCRIPT_DIR / "expression/fit_between_gene_expression.py",
        helper = EXPRESSION_SCRIPT_DIR / "expression/expression_model_helpers.py",
        data = R2C2_DATA / "r2c2/expression_model_data.tsv"
    output:
        coefficients = FUNCTIONAL_DIR / "between_gene/model_coefficients.tsv",
        support = FUNCTIONAL_DIR / "between_gene/model_support.tsv",
        summary = FUNCTIONAL_DIR / "between_gene/summary.txt"
    params: outdir = FUNCTIONAL_DIR / "between_gene"
    log: EXPRESSION_LOG_DIR / "fit_between_gene_expression.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --input {input.data:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """
