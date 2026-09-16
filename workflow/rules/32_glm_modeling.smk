# RNA-seq and R2C2 analyses promoted from the analysis inventory.
# Model runs and BAM scans are submitted through Snakemake by the user.
GLM_DIR = EXPRESSION_DIR / "glm_modeling"
FUNCTIONAL_DIR = EXPRESSION_DIR / "functional"
EXPRESSION_SCRIPT_DIR = PROJECT_ROOT / "scripts" / "expression"
R2C2_QUANT_DIR = Path(config["r2c2_quant_dir"])
RAW_SQANTI = [PROJECT_ROOT / "data" / f"sqanti3_output_{label}" / f"{label}_isoforms_classification.filtered.txt" for label in ("834", "1614", "1749")]
R2C2_QUANTS = [R2C2_QUANT_DIR / f"09092025_{label}_Isoforms.filtered.clean.quant" for label in ("834", "1614", "1749")]
EXPRESSION_GTFS = expand(ANNOTATIONS_DIR / "{sample}.gtf", sample=RNA_SEQ_SAMPLES)
EXPRESSION_COUNTS = expand(COUNTS_DIR / "{replicate}.counts.txt", replicate=get_all_replicates())
SPLICE_BOUNDARIES = expand(INTRONER_SPLICE_DIR / "{sample}.per_locus.tsv", sample=RNA_SEQ_SAMPLES)
EXPRESSION_BAMS = expand(ALIGNMENT_DIR / "bams" / "{replicate}.sorted.bam", replicate=get_all_replicates())
EXPRESSION_BAIS = expand(ALIGNMENT_DIR / "bams" / "{replicate}.sorted.bam.bai", replicate=get_all_replicates())
EXPRESSION_JUNCTIONS = expand(REGTOOLS_JUNCTION_DIR / "{replicate}.junctions.bed", replicate=get_all_replicates())
SHORT_READ_DATA = FUNCTIONAL_DIR / "short_read_data"
R2C2_DATA = FUNCTIONAL_DIR / "r2c2_data"
GAIN_MODEL_DIR = GLM_DIR / "within_gene_gain"
ISOFORM_MODEL_DIR = GLM_DIR / "reference_independent_isoform_change"

rule calculate_gene_gc_content:
    input:
        script = EXPRESSION_SCRIPT_DIR / "calculate_gc_content.py",
        gtfs = EXPRESSION_GTFS,
        assemblies = expand(ASSEMBLIES_DIR / "{sample}.vg_paths.fa", sample=RNA_SEQ_SAMPLES)
    output:
        file0 = GLM_DIR / "gc_content_by_gene_strain.csv"
    params:
        outdir = GLM_DIR
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "calculate_gene_gc_content.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --samples CCMP1545 RCC1614 RCC1749 --gtfs {input.gtfs:q} --assemblies {input.assemblies:q} --output {output.file0:q} > {log:q} 2>&1
        """

rule prepare_short_read_expression:
    input:
        script = EXPRESSION_SCRIPT_DIR / "negative_binomial_expression_pilot/prepare_expression_pilot_data.py",
        counts = EXPRESSION_COUNTS,
        gtfs = EXPRESSION_GTFS,
        features = ISOFORM_TURNOVER_DIR / "introner_features_by_gene_strain.tsv",
        gc = GLM_DIR / "gc_content_by_gene_strain.csv"
    output:
        file0 = SHORT_READ_DATA / "expression_model_data.tsv",
        file1 = SHORT_READ_DATA / "normalization_factors.tsv",
        file2 = SHORT_READ_DATA / "data_summary.txt"
    params:
        outdir = SHORT_READ_DATA
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "prepare_short_read_expression.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --project-root {PROJECT_ROOT:q} --output-dir {params.outdir:q} > {log:q} 2>&1
        """

rule prepare_r2c2_expression:
    input:
        script = EXPRESSION_SCRIPT_DIR / "r2c2_expression/prepare_r2c2_expression.py",
        data = SHORT_READ_DATA / "expression_model_data.tsv",
        quants = R2C2_QUANTS,
        sqanti = RAW_SQANTI,
        gtfs = EXPRESSION_GTFS,
        helpers = [EXPRESSION_SCRIPT_DIR / "negative_binomial_expression_pilot" / name for name in ("fit_expression_pilot_models.py", "prepare_expression_pilot_data.py")]
    output:
        file0 = R2C2_DATA / "r2c2/expression_model_data.tsv",
        file1 = R2C2_DATA / "r2c2/normalization_factors.tsv",
        file2 = R2C2_DATA / "short_read_matched/expression_model_data.tsv",
        file3 = R2C2_DATA / "short_read_matched/normalization_factors.tsv",
        file4 = R2C2_DATA / "mapping_audit.tsv",
        file5 = R2C2_DATA / "isoform_gene_crosswalk.tsv",
        file6 = R2C2_DATA / "excluded_mating_type_gene_ids.tsv",
        file7 = R2C2_DATA / "provenance.json"
    params:
        outdir = R2C2_DATA
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "prepare_r2c2_expression.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --input {input.data:q} --quant-dir {R2C2_QUANT_DIR:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule fit_within_gene_gain:
    input:
        script = EXPRESSION_SCRIPT_DIR / "r2c2_expression/fit_within_gene_gain.py",
        data = R2C2_DATA / "r2c2/expression_model_data.tsv",
        helper = EXPRESSION_SCRIPT_DIR / "negative_binomial_expression_pilot/fit_expression_pilot_models.py",
        colors = PROJECT_ROOT / "scripts/figure_color_guide.py",
        guide = PROJECT_ROOT / "master_figure_color_guide.tsv"
    output:
        file0 = GAIN_MODEL_DIR / "summary.txt",
        file1 = GAIN_MODEL_DIR / "gain_effect.tsv",
        file2 = GAIN_MODEL_DIR / "model_coefficients.tsv",
        file3 = GAIN_MODEL_DIR / "model_data.tsv",
        file4 = GAIN_MODEL_DIR / "model_support.tsv",
        file5 = GAIN_MODEL_DIR / "within_gene_gain.png"
    params:
        outdir = GAIN_MODEL_DIR
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "fit_within_gene_gain.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --input {input.data:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule fit_reference_independent_isoform_change:
    input:
        script = EXPRESSION_SCRIPT_DIR / "glm_modeling/fit_reference_independent_isoform_change.py",
        data = ISOFORM_TURNOVER_DIR / "isoform_introner_model_data.tsv"
    output:
        file0 = ISOFORM_MODEL_DIR / "summary.txt",
        file1 = ISOFORM_MODEL_DIR / "conditional_poisson_coefficients.tsv",
        file2 = ISOFORM_MODEL_DIR / "support_threshold_sensitivity.tsv",
        file3 = ISOFORM_MODEL_DIR / "within_between_negative_binomial.tsv",
        file4 = ISOFORM_MODEL_DIR / "within_between_negative_binomial_sensitivity.tsv",
        file5 = ISOFORM_MODEL_DIR / "pairwise_gene_changes.tsv",
        file6 = ISOFORM_MODEL_DIR / "pairwise_change_model.tsv",
        file7 = ISOFORM_MODEL_DIR / "pairwise_change_summary.tsv",
        file8 = ISOFORM_MODEL_DIR / "pairwise_isoform_change.png"
    params:
        outdir = ISOFORM_MODEL_DIR
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "fit_reference_independent_isoform_change.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --input {input.data:q} --output-dir {params.outdir:q} --minimum-long-read-count 10 > {log:q} 2>&1
        """

rule glm_modeling_complete:
    input:
        GAIN_MODEL_DIR / "summary.txt",
        GAIN_MODEL_DIR / "gain_effect.tsv",
        GAIN_MODEL_DIR / "model_coefficients.tsv",
        GAIN_MODEL_DIR / "model_data.tsv",
        GAIN_MODEL_DIR / "model_support.tsv",
        GAIN_MODEL_DIR / "within_gene_gain.png",
        ISOFORM_MODEL_DIR / "summary.txt",
        ISOFORM_MODEL_DIR / "conditional_poisson_coefficients.tsv",
        ISOFORM_MODEL_DIR / "support_threshold_sensitivity.tsv",
        ISOFORM_MODEL_DIR / "within_between_negative_binomial.tsv",
        ISOFORM_MODEL_DIR / "within_between_negative_binomial_sensitivity.tsv",
        ISOFORM_MODEL_DIR / "pairwise_gene_changes.tsv",
        ISOFORM_MODEL_DIR / "pairwise_change_model.tsv",
        ISOFORM_MODEL_DIR / "pairwise_change_summary.tsv",
        ISOFORM_MODEL_DIR / "pairwise_isoform_change.png"

rule poisson_glm_only:
    input:
        ISOFORM_MODEL_DIR / "summary.txt",
        ISOFORM_MODEL_DIR / "conditional_poisson_coefficients.tsv",
        ISOFORM_MODEL_DIR / "support_threshold_sensitivity.tsv",
        ISOFORM_MODEL_DIR / "within_between_negative_binomial.tsv",
        ISOFORM_MODEL_DIR / "within_between_negative_binomial_sensitivity.tsv",
        ISOFORM_MODEL_DIR / "pairwise_gene_changes.tsv",
        ISOFORM_MODEL_DIR / "pairwise_change_model.tsv",
        ISOFORM_MODEL_DIR / "pairwise_change_summary.tsv",
        ISOFORM_MODEL_DIR / "pairwise_isoform_change.png"

rule isoform_turnover_models:
    input:
        ISOFORM_MODEL_DIR / "summary.txt",
        ISOFORM_MODEL_DIR / "conditional_poisson_coefficients.tsv",
        ISOFORM_MODEL_DIR / "support_threshold_sensitivity.tsv",
        ISOFORM_MODEL_DIR / "within_between_negative_binomial.tsv",
        ISOFORM_MODEL_DIR / "within_between_negative_binomial_sensitivity.tsv",
        ISOFORM_MODEL_DIR / "pairwise_gene_changes.tsv",
        ISOFORM_MODEL_DIR / "pairwise_change_model.tsv",
        ISOFORM_MODEL_DIR / "pairwise_change_summary.tsv",
        ISOFORM_MODEL_DIR / "pairwise_isoform_change.png"

rule expression_glm_only:
    input:
        GAIN_MODEL_DIR / "summary.txt",
        GAIN_MODEL_DIR / "gain_effect.tsv",
        GAIN_MODEL_DIR / "model_coefficients.tsv",
        GAIN_MODEL_DIR / "model_data.tsv",
        GAIN_MODEL_DIR / "model_support.tsv",
        GAIN_MODEL_DIR / "within_gene_gain.png"

rule expression_models:
    input:
        GAIN_MODEL_DIR / "summary.txt",
        GAIN_MODEL_DIR / "gain_effect.tsv",
        GAIN_MODEL_DIR / "model_coefficients.tsv",
        GAIN_MODEL_DIR / "model_data.tsv",
        GAIN_MODEL_DIR / "model_support.tsv",
        GAIN_MODEL_DIR / "within_gene_gain.png",
        ISOFORM_MODEL_DIR / "summary.txt",
        ISOFORM_MODEL_DIR / "conditional_poisson_coefficients.tsv",
        ISOFORM_MODEL_DIR / "support_threshold_sensitivity.tsv",
        ISOFORM_MODEL_DIR / "within_between_negative_binomial.tsv",
        ISOFORM_MODEL_DIR / "within_between_negative_binomial_sensitivity.tsv",
        ISOFORM_MODEL_DIR / "pairwise_gene_changes.tsv",
        ISOFORM_MODEL_DIR / "pairwise_change_model.tsv",
        ISOFORM_MODEL_DIR / "pairwise_change_summary.tsv",
        ISOFORM_MODEL_DIR / "pairwise_isoform_change.png"
