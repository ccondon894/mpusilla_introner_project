# Candidate selection and transcript-oriented representative isoform figures (17, 20).
R2C2_ISOFORM_DIR = Path(config["r2c2_isoform_dir"])

rule select_isoform_candidates:
    input:
        script = EXPRESSION_SCRIPT_DIR / "isoform_analysis/find_pairwise_isoform_candidate_genes.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        sqanti = RAW_SQANTI[:2],
        isoform_gtfs = [R2C2_ISOFORM_DIR / name for name in ("09092025_834_Isoforms.filtered.clean.scaffold_corrected.sorted.gtf", "09092025_1614_Isoforms.filtered.clean.sorted.gtf")]
    output:
        file0 = FUNCTIONAL_DIR / "isoform_candidates" / "pairwise_1614_isoform_gene_candidates.tsv",
        file1 = FUNCTIONAL_DIR / "isoform_candidates" / "current_genotype_introner_vs_annotated_intron_size.tsv"
    params:
        outdir = FUNCTIONAL_DIR / "isoform_candidates"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "select_isoform_candidates.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --project-root {PROJECT_ROOT:q} --isoform-dir {R2C2_ISOFORM_DIR:q} --output {output.file0:q} --sidecar {output.file1:q} > {log:q} 2>&1
        """

rule plot_paper_isoform_candidate_1:
    input:
        script = EXPRESSION_SCRIPT_DIR / "isoform_analysis/rehaul_introner_isoform_figure.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        map = ISOFORM_TURNOVER_DIR / "locus_gene_map.tsv",
        gtfs = EXPRESSION_GTFS,
        nmd = SQANTI_DIR / "parsed_sqanti3_data.tsv",
        raw_sqanti = RAW_SQANTI,
        isoform_gtfs = [R2C2_ISOFORM_DIR / name for name in ("09092025_834_Isoforms.filtered.clean.scaffold_corrected.sorted.gtf", "09092025_1749_Isoforms.filtered.clean.sorted.gtf")],
        base = EXPRESSION_SCRIPT_DIR / "isoform_analysis/plot_introner_isoform_tracks.py",
        colors = PROJECT_ROOT / "master_figure_color_guide.tsv"
    output:
        file0 = FIGURES_DIR / "paper_isoform_tracks" / "isoform_candidate_1.png",
        file1 = FIGURES_DIR / "paper_isoform_tracks" / "isoform_candidate_1.summary.tsv"
    params:
        outdir = FIGURES_DIR / "paper_isoform_tracks"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "plot_paper_isoform_candidate_1.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --genotype-matrix {input.matrix:q} --locus-gene-map {input.map:q} --annotation-dir {ANNOTATIONS_DIR:q} --nmd-table {input.nmd:q} --isoform-dir {R2C2_ISOFORM_DIR:q} --sqanti-dir {PROJECT_ROOT:q}/data --gene MicpuC2.EuGene.0000010132.3.0.228 --single-output {output.file0:q} --x-padding-bp 100 --hide-track-separators --plain-introner-labels --collapse-matching-structures > {log:q} 2>&1
        """

rule plot_paper_isoform_candidate_2:
    input:
        script = EXPRESSION_SCRIPT_DIR / "isoform_analysis/rehaul_introner_isoform_figure.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        map = ISOFORM_TURNOVER_DIR / "locus_gene_map.tsv",
        gtfs = EXPRESSION_GTFS,
        nmd = SQANTI_DIR / "parsed_sqanti3_data.tsv",
        raw_sqanti = RAW_SQANTI,
        isoform_gtfs = [R2C2_ISOFORM_DIR / name for name in ("09092025_834_Isoforms.filtered.clean.scaffold_corrected.sorted.gtf", "09092025_1749_Isoforms.filtered.clean.sorted.gtf")],
        base = EXPRESSION_SCRIPT_DIR / "isoform_analysis/plot_introner_isoform_tracks.py",
        colors = PROJECT_ROOT / "master_figure_color_guide.tsv"
    output:
        file0 = FIGURES_DIR / "paper_isoform_tracks" / "isoform_candidate_2.png",
        file1 = FIGURES_DIR / "paper_isoform_tracks" / "isoform_candidate_2.summary.tsv"
    params:
        outdir = FIGURES_DIR / "paper_isoform_tracks"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "plot_paper_isoform_candidate_2.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --genotype-matrix {input.matrix:q} --locus-gene-map {input.map:q} --annotation-dir {ANNOTATIONS_DIR:q} --nmd-table {input.nmd:q} --isoform-dir {R2C2_ISOFORM_DIR:q} --sqanti-dir {PROJECT_ROOT:q}/data --gene MicpuC2.estExt_Genewise1.C_100540.3.0.228 --single-output {output.file0:q} --x-padding-bp 100 --hide-track-separators --plain-introner-labels > {log:q} 2>&1
        """

rule paper_isoform_figures:
    input:
        FIGURES_DIR / "paper_isoform_tracks" / "isoform_candidate_1.png",
        FIGURES_DIR / "paper_isoform_tracks" / "isoform_candidate_1.summary.tsv",
        FIGURES_DIR / "paper_isoform_tracks" / "isoform_candidate_2.png",
        FIGURES_DIR / "paper_isoform_tracks" / "isoform_candidate_2.summary.tsv"

rule plot_functional_associations:
    input:
        script = EXPRESSION_SCRIPT_DIR / "functional_associations_figure/plot_functional_associations_data_forward_r2c2.py",
        helper = EXPRESSION_SCRIPT_DIR / "functional_associations_figure/plot_functional_associations_variants.py",
        colors = PROJECT_ROOT / "master_figure_color_guide.tsv",
        isoform = ISOFORM_MODEL_DIR / "conditional_poisson_coefficients.tsv",
        robustness = ISOFORM_MODEL_DIR / "within_between_negative_binomial.tsv",
        between = FUNCTIONAL_DIR / "r2c2_between_gene/model_coefficients.tsv",
        gain = GAIN_MODEL_DIR / "gain_effect.tsv",
        nmd = FUNCTIONAL_DIR / "adjusted_nmd/model_coefficients.tsv",
        nmd_data = FUNCTIONAL_DIR / "adjusted_nmd/nmd_gene_strain_data.tsv",
        splicing = FUNCTIONAL_DIR / "corrected_retention_psi/per_locus_scores.tsv",
        tests = FUNCTIONAL_DIR / "corrected_retention_psi/evidence_guarded_class_tests.tsv"
    output:
        figure = FIGURES_DIR / "functional_associations.png",
        statistics = FIGURES_DIR / "functional_associations.panel_b.tsv"
    log: EXPRESSION_LOG_DIR / "plot_functional_associations.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p $(dirname {output.figure:q}) $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --project-root {PROJECT_ROOT:q} --output {output.figure:q} > {log:q} 2>&1
        """

rule audit_r2c2_quantification:
    input:
        script = EXPRESSION_SCRIPT_DIR / "r2c2_expression/compare_834_quant_versions.py",
        source = R2C2_QUANT_DIR / "09092025_834_Isoforms.filtered.clean.quant",
        candidates = sorted(R2C2_ISOFORM_DIR.glob("*834*clean.quant*")) + sorted(R2C2_ISOFORM_DIR.glob("new_834*/**/Isoforms.filtered.clean.quant")),
        classifications = sorted((R2C2_ISOFORM_DIR / "sqanti3_output_834").glob("834_isoforms_classification.filtered.txt*")),
        sqanti = RAW_SQANTI[0],
        crosswalk = R2C2_DATA / "isoform_gene_crosswalk.tsv"
    output:
        summary = FUNCTIONAL_DIR / "quantification_audit/summary.txt",
        files = FUNCTIONAL_DIR / "quantification_audit/file_comparison.tsv",
        assignments = FUNCTIONAL_DIR / "quantification_audit/all_isoform_assignments.tsv",
        changes = FUNCTIONAL_DIR / "quantification_audit/changed_gene_labels.tsv",
        sqanti = FUNCTIONAL_DIR / "quantification_audit/sqanti_comparison.tsv"
    params:
        outdir = FUNCTIONAL_DIR / "quantification_audit"
    log:
        EXPRESSION_LOG_DIR / "audit_r2c2_quantification.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        python {input.script:q} --project-root {PROJECT_ROOT:q} \
            --comparison-dir {R2C2_ISOFORM_DIR:q} --source {input.source:q} \
            --crosswalk {input.crosswalk:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule expression_analysis_complete:
    input:
        rules.expression_models.input,
        FUNCTIONAL_DIR / "introner_length_classes" / "length_tests.tsv",
        FUNCTIONAL_DIR / "introner_splice_site_classes" / "annotation_class_tests.tsv",
        FUNCTIONAL_DIR / "introner_splicing_classes" / "statistical_tests.tsv",
        FUNCTIONAL_DIR / "rare_introner_splicing_classes" / "statistical_tests.tsv",
        FUNCTIONAL_DIR / "rnaseq_expression_classes" / "continuous_expression_tests.tsv",
        FUNCTIONAL_DIR / "isoform_qualified_expression_classes" / "statistical_tests.tsv",
        FUNCTIONAL_DIR / "rare_introner_expression_classes" / "statistical_tests.tsv",
        FUNCTIONAL_DIR / "conventional_intron_splicing" / "per_locus_scores.tsv",
        FUNCTIONAL_DIR / "corrected_retention_psi" / "per_locus_scores.tsv",
        FUNCTIONAL_DIR / "corrected_retention_psi" / "evidence_guarded_class_tests.tsv",
        FUNCTIONAL_DIR / "adjusted_retention_rates" / "adjusted_class_retention_rates.tsv",
        FUNCTIONAL_DIR / "depth_retention_psi" / "per_locus_depth_psi.tsv",
        FUNCTIONAL_DIR / "adjusted_nmd" / "model_coefficients.tsv",
        FUNCTIONAL_DIR / "group1_polymorphic_expression_pilot/models" / "test_a_coefficients.tsv",
        FUNCTIONAL_DIR / "r2c2_total_burden_expression_pilot/models" / "test_a_coefficients.tsv",
        FUNCTIONAL_DIR / "short_read_models" / "model_coefficients.tsv",
        FUNCTIONAL_DIR / "paired_introner_expression" / "model_coefficients.tsv",
        FUNCTIONAL_DIR / "current_burden_r2c2" / "model_coefficients.tsv",
        FUNCTIONAL_DIR / "current_burden_short_read_matched" / "model_coefficients.tsv",
        FUNCTIONAL_DIR / "r2c2_group1_gain" / "model_coefficients.tsv",
        FUNCTIONAL_DIR / "single_locus_expression" / "model_coefficients.tsv",
        rules.paper_isoform_figures.input,
        FUNCTIONAL_DIR / "isoform_candidates/pairwise_1614_isoform_gene_candidates.tsv",
        FIGURES_DIR / "functional_associations.png",
        FUNCTIONAL_DIR / "r2c2_between_gene/model_coefficients.tsv"
