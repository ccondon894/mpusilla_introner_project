# Additional RNA-seq/R2C2 occupancy and paired-expression analyses.

rule prepare_group1_polymorphic_expression_pilot:
    input:
        script = EXPRESSION_SCRIPT_DIR / "group1_polymorphic_expression_pilot/prepare_group1_expression_tables.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        mapping = ISOFORM_TURNOVER_DIR / "locus_gene_map.tsv",
        gtfs = EXPRESSION_GTFS,
        gc = GLM_DIR / "gc_content_by_gene_strain.csv",
        helper = EXPRESSION_SCRIPT_DIR / "negative_binomial_expression_pilot/prepare_expression_pilot_data.py",
        counts = EXPRESSION_COUNTS
    output:
        file0 = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/data" / "expression_gene_strain.tsv",
        file1 = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/data" / "gene_pair_table.tsv",
        file2 = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/data" / "locus_status.tsv",
        file3 = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/data" / "data_summary.txt"
    params:
        outdir = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/data"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "prepare_group1_polymorphic_expression_pilot.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --project-root {PROJECT_ROOT:q} --output-dir {params.outdir:q} > {log:q} 2>&1
        """

rule fit_group1_polymorphic_expression_pilot:
    input:
        script = EXPRESSION_SCRIPT_DIR / "group1_polymorphic_expression_pilot/fit_group1_expression_models.py",
        expression = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/data" / "expression_gene_strain.tsv",
        pairs = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/data" / "gene_pair_table.tsv",
        loci = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/data" / "locus_status.tsv",
        colors = PROJECT_ROOT / "scripts/figure_color_guide.py",
        guide = PROJECT_ROOT / "master_figure_color_guide.tsv"
    output:
        file0 = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/models" / "test_a_coefficients.tsv",
        file1 = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/models" / "test_b_coefficients.tsv",
        file2 = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/models" / "model_summary.txt",
        file3 = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/models" / "validation_summary.txt"
    params:
        outdir = FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/models"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "fit_group1_polymorphic_expression_pilot.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --data-dir {params.outdir:q}/../data --output-dir {params.outdir:q} > {log:q} 2>&1
        """

rule prepare_r2c2_total_burden_expression_pilot:
    input:
        script = EXPRESSION_SCRIPT_DIR / "r2c2_total_burden_expression_pilot/prepare_r2c2_burden_tables.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        mapping = ISOFORM_TURNOVER_DIR / "locus_gene_map.tsv",
        gtfs = EXPRESSION_GTFS,
        gc = GLM_DIR / "gc_content_by_gene_strain.csv",
        helper = EXPRESSION_SCRIPT_DIR / "negative_binomial_expression_pilot/prepare_expression_pilot_data.py",
        quants = R2C2_QUANTS,
        sqanti = RAW_SQANTI
    output:
        file0 = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/data" / "expression_gene_strain.tsv",
        file1 = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/data" / "gene_pair_table.tsv",
        file2 = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/data" / "locus_status.tsv",
        file3 = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/data" / "data_summary.txt"
    params:
        outdir = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/data"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "prepare_r2c2_total_burden_expression_pilot.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --project-root {PROJECT_ROOT:q} --output-dir {params.outdir:q} --quant-dir {R2C2_QUANT_DIR:q} > {log:q} 2>&1
        """

rule fit_r2c2_total_burden_expression_pilot:
    input:
        script = EXPRESSION_SCRIPT_DIR / "r2c2_total_burden_expression_pilot/fit_r2c2_burden_models.py",
        expression = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/data" / "expression_gene_strain.tsv",
        pairs = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/data" / "gene_pair_table.tsv",
        loci = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/data" / "locus_status.tsv",
        colors = PROJECT_ROOT / "scripts/figure_color_guide.py",
        guide = PROJECT_ROOT / "master_figure_color_guide.tsv"
    output:
        file0 = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/models" / "test_a_coefficients.tsv",
        file1 = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/models" / "test_b_coefficients.tsv",
        file2 = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/models" / "model_summary.txt",
        file3 = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/models" / "validation_summary.txt"
    params:
        outdir = FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/models"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "fit_r2c2_total_burden_expression_pilot.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --data-dir {params.outdir:q}/../data --output-dir {params.outdir:q} > {log:q} 2>&1
        """

rule short_read_expression_models:
    input:
        script = EXPRESSION_SCRIPT_DIR / "negative_binomial_expression_pilot/fit_expression_pilot_models.py",
        data = SHORT_READ_DATA / "expression_model_data.tsv",
        normalization = SHORT_READ_DATA / "normalization_factors.tsv"
    output:
        file0 = FUNCTIONAL_DIR / "short_read_models" / "model_coefficients.tsv",
        file1 = FUNCTIONAL_DIR / "short_read_models" / "model_summary.txt",
        file2 = FUNCTIONAL_DIR / "short_read_models" / "model_support.tsv"
    params:
        outdir = FUNCTIONAL_DIR / "short_read_models"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "short_read_expression_models.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --input {input.data:q} --normalization {input.normalization:q} --output-dir {params.outdir:q} > {log:q} 2>&1
        """

rule paired_introner_expression:
    input:
        script = EXPRESSION_SCRIPT_DIR / "paired_introner_expression/analyze_paired_expression.py",
        data = SHORT_READ_DATA / "expression_model_data.tsv",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        mapping = ISOFORM_TURNOVER_DIR / "locus_gene_map.tsv",
        events = ISOFORM_TURNOVER_DIR / "locus_level_reference_comparisons.tsv",
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        helpers = [EXPRESSION_SCRIPT_DIR / "negative_binomial_expression_pilot" / name for name in ("prepare_expression_pilot_data.py", "fit_expression_pilot_models.py")]
    output:
        file0 = FUNCTIONAL_DIR / "paired_introner_expression" / "model_coefficients.tsv",
        file1 = FUNCTIONAL_DIR / "paired_introner_expression" / "validation_summary.txt",
        file2 = FUNCTIONAL_DIR / "paired_introner_expression" / "results_summary.txt",
        file3 = FUNCTIONAL_DIR / "paired_introner_expression" / "paired_expression_effects.png"
    params:
        outdir = FUNCTIONAL_DIR / "paired_introner_expression"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "paired_introner_expression.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --project-root {PROJECT_ROOT:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule current_burden_r2c2:
    input:
        script = EXPRESSION_SCRIPT_DIR / "r2c2_expression/run_current_burden_models.py",
        data = R2C2_DATA / "r2c2/expression_model_data.tsv",
        helpers = [EXPRESSION_SCRIPT_DIR / "r2c2_expression/prepare_r2c2_expression.py", EXPRESSION_SCRIPT_DIR / "negative_binomial_expression_pilot/fit_expression_pilot_models.py"]
    output:
        file0 = FUNCTIONAL_DIR / "current_burden_r2c2" / "model_coefficients.tsv",
        file1 = FUNCTIONAL_DIR / "current_burden_r2c2" / "burden_effects.tsv",
        file2 = FUNCTIONAL_DIR / "current_burden_r2c2" / "model_support.tsv",
        file3 = FUNCTIONAL_DIR / "current_burden_r2c2" / "summary.txt",
        file4 = FUNCTIONAL_DIR / "current_burden_r2c2" / "validation_summary.txt"
    params:
        outdir = FUNCTIONAL_DIR / "current_burden_r2c2"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "current_burden_r2c2.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --data-dir {R2C2_DATA:q} --count-source r2c2 --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule current_burden_short_read_matched:
    input:
        script = EXPRESSION_SCRIPT_DIR / "r2c2_expression/run_current_burden_models.py",
        data = R2C2_DATA / "short_read_matched/expression_model_data.tsv",
        helpers = [EXPRESSION_SCRIPT_DIR / "r2c2_expression/prepare_r2c2_expression.py", EXPRESSION_SCRIPT_DIR / "negative_binomial_expression_pilot/fit_expression_pilot_models.py"]
    output:
        file0 = FUNCTIONAL_DIR / "current_burden_short_read_matched" / "model_coefficients.tsv",
        file1 = FUNCTIONAL_DIR / "current_burden_short_read_matched" / "burden_effects.tsv",
        file2 = FUNCTIONAL_DIR / "current_burden_short_read_matched" / "model_support.tsv",
        file3 = FUNCTIONAL_DIR / "current_burden_short_read_matched" / "summary.txt",
        file4 = FUNCTIONAL_DIR / "current_burden_short_read_matched" / "validation_summary.txt"
    params:
        outdir = FUNCTIONAL_DIR / "current_burden_short_read_matched"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "current_burden_short_read_matched.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --data-dir {R2C2_DATA:q} --count-source short_read_matched --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule r2c2_group1_gain:
    input:
        script = EXPRESSION_SCRIPT_DIR / "r2c2_expression/run_group1_gain_model.py",
        data = R2C2_DATA / "r2c2/expression_model_data.tsv",
        coefficients = GAIN_MODEL_DIR / "model_coefficients.tsv",
        helpers = [EXPRESSION_SCRIPT_DIR / "r2c2_expression/prepare_r2c2_expression.py", EXPRESSION_SCRIPT_DIR / "negative_binomial_expression_pilot/fit_expression_pilot_models.py"]
    output:
        file0 = FUNCTIONAL_DIR / "r2c2_group1_gain" / "model_coefficients.tsv",
        file1 = FUNCTIONAL_DIR / "r2c2_group1_gain" / "gain_loss_tests.tsv",
        file2 = FUNCTIONAL_DIR / "r2c2_group1_gain" / "model_support.tsv",
        file3 = FUNCTIONAL_DIR / "r2c2_group1_gain" / "summary.txt"
    params:
        outdir = FUNCTIONAL_DIR / "r2c2_group1_gain"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "r2c2_group1_gain.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --data-dir {R2C2_DATA:q} --gain-coefficients {input.coefficients:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule single_locus_expression:
    input:
        script = EXPRESSION_SCRIPT_DIR / "r2c2_expression/run_single_locus_three_strains.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        events = ISOFORM_TURNOVER_DIR / "locus_level_reference_comparisons.tsv",
        map = ISOFORM_TURNOVER_DIR / "locus_gene_map.tsv",
        long_data = R2C2_DATA / "r2c2/expression_model_data.tsv",
        short_data = R2C2_DATA / "short_read_matched/expression_model_data.tsv",
        mt = R2C2_DATA / "excluded_mating_type_gene_ids.tsv",
        helpers = [EXPRESSION_SCRIPT_DIR / "r2c2_expression/prepare_r2c2_expression.py", EXPRESSION_SCRIPT_DIR / "negative_binomial_expression_pilot/fit_expression_pilot_models.py"]
    output:
        file0 = FUNCTIONAL_DIR / "single_locus_expression" / "model_coefficients.tsv",
        file1 = FUNCTIONAL_DIR / "single_locus_expression" / "presence_effects.tsv",
        file2 = FUNCTIONAL_DIR / "single_locus_expression" / "model_support.tsv",
        file3 = FUNCTIONAL_DIR / "single_locus_expression" / "summary.txt"
    params:
        outdir = FUNCTIONAL_DIR / "single_locus_expression"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "single_locus_expression.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --data-dir {R2C2_DATA:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """
