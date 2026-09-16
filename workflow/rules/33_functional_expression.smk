# Expression classes, splicing, retention and NMD (inventory entries 5–13, 25).

rule introner_length_classes:
    input:
        script = EXPRESSION_SCRIPT_DIR / "functional/compare_fixed_polymorphic_introner_lengths.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        boundaries = SPLICE_BOUNDARIES
    output:
        file0 = FUNCTIONAL_DIR / "introner_length_classes" / "length_tests.tsv",
        file1 = FUNCTIONAL_DIR / "introner_length_classes" / "length_summary.tsv",
        file2 = FUNCTIONAL_DIR / "introner_length_classes" / "introner_length_comparison.png",
        file3 = FUNCTIONAL_DIR / "introner_length_classes" / "RESULTS.md"
    params:
        outdir = FUNCTIONAL_DIR / "introner_length_classes"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "introner_length_classes.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --matrix {input.matrix:q} --boundary-dir {INTRONER_SPLICE_DIR:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule introner_splice_site_classes:
    input:
        script = EXPRESSION_SCRIPT_DIR / "functional/compare_fixed_polymorphic_splice_sites.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        boundaries = SPLICE_BOUNDARIES,
        assemblies = expand(ASSEMBLIES_DIR / "{sample}.vg_paths.fa", sample=RNA_SEQ_SAMPLES)
    output:
        file0 = FUNCTIONAL_DIR / "introner_splice_site_classes" / "annotation_class_tests.tsv",
        file1 = FUNCTIONAL_DIR / "introner_splice_site_classes" / "rna_supported_motif_tests.tsv",
        file2 = FUNCTIONAL_DIR / "introner_splice_site_classes" / "splice_site_class_comparison.png",
        file3 = FUNCTIONAL_DIR / "introner_splice_site_classes" / "RESULTS.md"
    params:
        outdir = FUNCTIONAL_DIR / "introner_splice_site_classes"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "introner_splice_site_classes.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --matrix {input.matrix:q} --boundary-dir {INTRONER_SPLICE_DIR:q} --outdir {params.outdir:q} --assembly-dir {ASSEMBLIES_DIR:q} > {log:q} 2>&1
        """

rule introner_splicing_classes:
    input:
        script = EXPRESSION_SCRIPT_DIR / "functional/compare_introner_splicing_efficiency.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        boundaries = SPLICE_BOUNDARIES
    output:
        file0 = FUNCTIONAL_DIR / "introner_splicing_classes" / "statistical_tests.tsv",
        file1 = FUNCTIONAL_DIR / "introner_splicing_classes" / "per_locus_classified.tsv",
        file2 = FUNCTIONAL_DIR / "introner_splicing_classes" / "splicing_efficiency_psi_distributions.png"
    params:
        outdir = FUNCTIONAL_DIR / "introner_splicing_classes"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "introner_splicing_classes.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --matrix {input.matrix:q} --psi-dir {INTRONER_SPLICE_DIR:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule rare_introner_splicing_classes:
    input:
        script = EXPRESSION_SCRIPT_DIR / "functional/compare_introner_splicing_efficiency.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        boundaries = SPLICE_BOUNDARIES
    output:
        file0 = FUNCTIONAL_DIR / "rare_introner_splicing_classes" / "statistical_tests.tsv",
        file1 = FUNCTIONAL_DIR / "rare_introner_splicing_classes" / "per_locus_classified.tsv",
        file2 = FUNCTIONAL_DIR / "rare_introner_splicing_classes" / "splicing_efficiency_psi_distributions.png"
    params:
        outdir = FUNCTIONAL_DIR / "rare_introner_splicing_classes"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "rare_introner_splicing_classes.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --matrix {input.matrix:q} --psi-dir {INTRONER_SPLICE_DIR:q} --outdir {params.outdir:q} --polymorphic-min-total-present 1 --polymorphic-max-total-present 3 > {log:q} 2>&1
        """

rule rnaseq_expression_classes:
    input:
        script = EXPRESSION_SCRIPT_DIR / "functional/compare_introner_expression_by_class_counts.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        counts = COUNTS_DIR / "merged_counts_matrix.csv"
    output:
        file0 = FUNCTIONAL_DIR / "rnaseq_expression_classes" / "continuous_expression_tests.tsv",
        file1 = FUNCTIONAL_DIR / "rnaseq_expression_classes" / "low_expression_tests.tsv",
        file2 = FUNCTIONAL_DIR / "rnaseq_expression_classes" / "negative_binomial_gee.tsv",
        file3 = FUNCTIONAL_DIR / "rnaseq_expression_classes" / "gene_expression_distributions.png",
        file4 = FUNCTIONAL_DIR / "rnaseq_expression_classes" / "RESULTS.md"
    params:
        outdir = FUNCTIONAL_DIR / "rnaseq_expression_classes"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "rnaseq_expression_classes.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --matrix {input.matrix:q} --counts {input.counts:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule prepare_isoform_qualified_expression:
    input:
        script = EXPRESSION_SCRIPT_DIR / "functional/prepare_isoform_qualified_expression.py",
        sqanti = SQANTI_DIR / "parsed_sqanti3_data.tsv",
        counts = COUNTS_DIR / "merged_counts_matrix.csv",
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf"
    output:
        file0 = FUNCTIONAL_DIR / "isoform_qualified_expression" / "expression.csv"
    params:
        outdir = FUNCTIONAL_DIR / "isoform_qualified_expression"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "prepare_isoform_qualified_expression.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --sqanti {input.sqanti:q} --counts {input.counts:q} --gtf {input.gtf:q} --output {output.file0:q} > {log:q} 2>&1
        """

rule isoform_qualified_expression_classes:
    input:
        script = EXPRESSION_SCRIPT_DIR / "functional/compare_introner_expression_by_class.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        expression = FUNCTIONAL_DIR / "isoform_qualified_expression/expression.csv"
    output:
        file0 = FUNCTIONAL_DIR / "isoform_qualified_expression_classes" / "statistical_tests.tsv",
        file1 = FUNCTIONAL_DIR / "isoform_qualified_expression_classes" / "gene_class_expression.tsv",
        file2 = FUNCTIONAL_DIR / "isoform_qualified_expression_classes" / "gene_expression_distributions.png"
    params:
        outdir = FUNCTIONAL_DIR / "isoform_qualified_expression_classes"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "isoform_qualified_expression_classes.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --matrix {input.matrix:q} --expression {input.expression:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule rare_introner_expression_classes:
    input:
        script = EXPRESSION_SCRIPT_DIR / "functional/compare_introner_expression_by_class.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        expression = FUNCTIONAL_DIR / "isoform_qualified_expression/expression.csv"
    output:
        file0 = FUNCTIONAL_DIR / "rare_introner_expression_classes" / "statistical_tests.tsv",
        file1 = FUNCTIONAL_DIR / "rare_introner_expression_classes" / "gene_class_expression.tsv",
        file2 = FUNCTIONAL_DIR / "rare_introner_expression_classes" / "gene_expression_distributions.png"
    params:
        outdir = FUNCTIONAL_DIR / "rare_introner_expression_classes"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "rare_introner_expression_classes.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --matrix {input.matrix:q} --expression {input.expression:q} --outdir {params.outdir:q} --polymorphic-min-total-present 1 --polymorphic-max-total-present 3 > {log:q} 2>&1
        """

rule conventional_intron_splicing:
    input:
        script = EXPRESSION_SCRIPT_DIR / "introner_vs_non_introner_splicing_test/compare_to_conventional_introns.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        gtfs = EXPRESSION_GTFS,
        junctions = EXPRESSION_JUNCTIONS,
        boundaries = SPLICE_BOUNDARIES,
        bams = EXPRESSION_BAMS,
        bais = EXPRESSION_BAIS,
        helper = EXPRESSION_SCRIPT_DIR / "introner_vs_non_introner_splicing_test/compare_to_conventional_introns.py"
    output:
        file0 = FUNCTIONAL_DIR / "conventional_intron_splicing" / "per_locus_scores.tsv",
        file1 = FUNCTIONAL_DIR / "conventional_intron_splicing" / "gene_cluster_aware_tests.tsv",
        file2 = FUNCTIONAL_DIR / "conventional_intron_splicing" / "RESULTS.md"
    params:
        outdir = FUNCTIONAL_DIR / "conventional_intron_splicing"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "conventional_intron_splicing.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --matrix {input.matrix:q} --annotations {ANNOTATIONS_DIR:q} --splice-dir {SPLICE_JUNCTION_DIR:q} --bam-dir {ALIGNMENT_DIR:q}/bams --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule corrected_retention_psi:
    input:
        script = EXPRESSION_SCRIPT_DIR / "introner_vs_non_introner_splicing_test/corrected_retention_psi.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        gtfs = EXPRESSION_GTFS,
        junctions = EXPRESSION_JUNCTIONS,
        boundaries = SPLICE_BOUNDARIES,
        bams = EXPRESSION_BAMS,
        bais = EXPRESSION_BAIS,
        helper = EXPRESSION_SCRIPT_DIR / "introner_vs_non_introner_splicing_test/compare_to_conventional_introns.py"
    output:
        file0 = FUNCTIONAL_DIR / "corrected_retention_psi" / "per_locus_scores.tsv",
        file1 = FUNCTIONAL_DIR / "corrected_retention_psi" / "per_locus_replicate_U5_U3_S.tsv",
        file2 = FUNCTIONAL_DIR / "corrected_retention_psi" / "boundary_resolver_audit.tsv",
        file3 = FUNCTIONAL_DIR / "corrected_retention_psi" / "validation_checks.tsv"
    params:
        outdir = FUNCTIONAL_DIR / "corrected_retention_psi"
    threads: 4
    log:
        EXPRESSION_LOG_DIR / "corrected_retention_psi.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --matrix {input.matrix:q} --annotations {ANNOTATIONS_DIR:q} --splice-dir {SPLICE_JUNCTION_DIR:q} --bam-dir {ALIGNMENT_DIR:q}/bams --outdir {params.outdir:q} --workers {threads} > {log:q} 2>&1
        """

rule summarize_corrected_retention_psi:
    input:
        script = EXPRESSION_SCRIPT_DIR / "introner_vs_non_introner_splicing_test/summarize_corrected_retention_psi.py",
        loci = FUNCTIONAL_DIR / "corrected_retention_psi/per_locus_scores.tsv",
        replicates = FUNCTIONAL_DIR / "corrected_retention_psi/per_locus_replicate_U5_U3_S.tsv",
        resolver = FUNCTIONAL_DIR / "corrected_retention_psi/boundary_resolver_audit.tsv"
    output:
        file0 = FUNCTIONAL_DIR / "corrected_retention_psi" / "evidence_guarded_class_tests.tsv",
        file1 = FUNCTIONAL_DIR / "corrected_retention_psi" / "PRIMARY_RESULTS.md"
    params:
        outdir = FUNCTIONAL_DIR / "corrected_retention_psi"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "summarize_corrected_retention_psi.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule adjusted_retention_rates:
    input:
        script = EXPRESSION_SCRIPT_DIR / "functional/fit_adjusted_retention_rates.py",
        counts = COUNTS_DIR / "merged_counts_matrix.csv",
        psi = FUNCTIONAL_DIR / "corrected_retention_psi/per_locus_replicate_U5_U3_S.tsv"
    output:
        file0 = FUNCTIONAL_DIR / "adjusted_retention_rates" / "adjusted_class_retention_rates.tsv",
        file1 = FUNCTIONAL_DIR / "adjusted_retention_rates" / "adjusted_class_retention_contrasts.tsv",
        file2 = FUNCTIONAL_DIR / "adjusted_retention_rates" / "model_summary.txt"
    params:
        outdir = FUNCTIONAL_DIR / "adjusted_retention_rates"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "adjusted_retention_rates.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --counts {input.counts:q} --psi-counts {input.psi:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

rule depth_retention_psi:
    input:
        script = EXPRESSION_SCRIPT_DIR / "introner_vs_non_introner_splicing_test/paper_depth_psi.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        gtfs = EXPRESSION_GTFS,
        junctions = EXPRESSION_JUNCTIONS,
        boundaries = SPLICE_BOUNDARIES,
        bams = EXPRESSION_BAMS,
        bais = EXPRESSION_BAIS,
        helper = EXPRESSION_SCRIPT_DIR / "introner_vs_non_introner_splicing_test/compare_to_conventional_introns.py",
        helper2 = EXPRESSION_SCRIPT_DIR / "introner_vs_non_introner_splicing_test/corrected_retention_psi.py"
    output:
        file0 = FUNCTIONAL_DIR / "depth_retention_psi" / "per_locus_depth_psi.tsv",
        file1 = FUNCTIONAL_DIR / "depth_retention_psi" / "RESULTS.md"
    params:
        outdir = FUNCTIONAL_DIR / "depth_retention_psi"
    threads: 4
    log:
        EXPRESSION_LOG_DIR / "depth_retention_psi.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --matrix {input.matrix:q} --annotations {ANNOTATIONS_DIR:q} --splice-dir {SPLICE_JUNCTION_DIR:q} --bam-dir {ALIGNMENT_DIR:q}/bams --outdir {params.outdir:q} --workers {threads} > {log:q} 2>&1
        """

rule adjusted_nmd_models:
    input:
        script = EXPRESSION_SCRIPT_DIR / "nmd_introner_adjusted/fit_adjusted_nmd_models.py",
        sqanti = RAW_SQANTI,
        gtfs = EXPRESSION_GTFS,
        architecture = ISOFORM_TURNOVER_DIR / "isoform_introner_model_data.tsv",
        map = ISOFORM_TURNOVER_DIR / "locus_gene_map.tsv",
        crosswalk = ISOFORM_TURNOVER_DIR / "sample_gene_crosswalk.tsv"
    output:
        file0 = FUNCTIONAL_DIR / "adjusted_nmd" / "model_coefficients.tsv",
        file1 = FUNCTIONAL_DIR / "adjusted_nmd" / "summary.txt",
        file2 = FUNCTIONAL_DIR / "adjusted_nmd" / "validation_summary.txt",
        data = FUNCTIONAL_DIR / "adjusted_nmd" / "nmd_gene_strain_data.tsv"
    params:
        outdir = FUNCTIONAL_DIR / "adjusted_nmd"
    threads: 1
    log:
        EXPRESSION_LOG_DIR / "adjusted_nmd_models.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp
        export MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --project-root {PROJECT_ROOT:q} --output-dir {params.outdir:q} > {log:q} 2>&1
        """
