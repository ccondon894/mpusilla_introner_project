# ============================================================
# 60_random_forest_classifier.smk - Random Forest Classifier
# ============================================================
#
# Present/absent introner locus classification and cross-strain
# pre-insertion scoring. Ports analysis/random_forest_classifier_project.
#
# Key analyses:
# - CCMP1545 within-sample leave-one-family-out
# - CCMP1545 train -> RCC1749 present/absent transfer
# - CCMP1545 train -> RCC1749 empty-ortholog pre-insertion scoring
# - Same-gene GC-matched control scoring and k-mer permutation
#
# ============================================================

from pathlib import Path
import shlex


# ============================================================
# CONFIGURATION
# ============================================================

RF_CONFIG = config.get("random_forest", {})
RF_TRAIN_SAMPLE = RF_CONFIG.get("train_sample", REFERENCE)
RF_TEST_SAMPLE = RF_CONFIG.get("test_sample", config["samples"]["reference"]["group2"])
RF_WINDOWS = RF_CONFIG.get("windows", [25, 50, 100])
RF_WINDOWS_STR = ",".join(str(window) for window in RF_WINDOWS)
RF_WINDOWS_TAG = "w" + "_".join(str(window) for window in RF_WINDOWS)
RF_MASK = RF_CONFIG.get("junction_mask_bp", 3)
RF_BOUNDARY_MAX_DELTA = RF_CONFIG.get("boundary_max_abs_delta", 25)
RF_EXCLUDED_CONTIGS = RF_CONFIG.get("excluded_contigs", {})
RF_ORIENTATION = RF_CONFIG.get("orientation", {})
RF_ENRICHMENT = RF_CONFIG.get("enrichment", {})
RF_MODELING = RF_CONFIG.get("modeling", {})
RF_FEATURE_SETS = RF_CONFIG.get("feature_sets", {})

RF_G1 = RF_ORIENTATION.get("mummer_pair", [RF_TRAIN_SAMPLE, RF_TEST_SAMPLE])[0]
RF_G2 = RF_ORIENTATION.get("mummer_pair", [RF_TRAIN_SAMPLE, RF_TEST_SAMPLE])[1]
RF_ORIENTATION_THRESHOLD = RF_ORIENTATION.get("threshold", 0.8)

RF_MICROC_RES = RF_ENRICHMENT.get("microc_resolution", 1000)
RF_MICROC_WINDOW = RF_ENRICHMENT.get("microc_window", 25000)
RF_MICROC_FLANK_WINDOWS = RF_ENRICHMENT.get("microc_flank_windows", [25, 50, 250, 500])
RF_MICROC_FLANK_STR = ",".join(str(window) for window in RF_MICROC_FLANK_WINDOWS)
RF_MICROC_PREFIX = f"{RF_MICROC_RES}bp_{RF_MICROC_WINDOW}bp"
RF_EXPR_PREFIX = RF_ENRICHMENT.get("expression_sample_prefix", "834")
RF_INCLUDE_EXPRESSION = RF_ENRICHMENT.get("include_expression", True)
RF_INCLUDE_PSI = RF_ENRICHMENT.get("include_psi", False)

RF_DIR = PROJECT_ROOT / config["output"].get("random_forest", "results/random_forest")
RF_FEATURES_DIR = RF_DIR / "features"
RF_MODELS_DIR = RF_DIR / "models"
RF_LOG_DIR = RF_DIR / "logs"
RF_SCRIPTS = PROJECT_ROOT / "scripts" / "random_forest"

MUMMER_COORDS = GENOME_ALIGNMENT_DIR / "mummer" / f"{RF_G1}_vs_{RF_G2}.coords"
BOUNDARY_SUPPORT_DIR = EXPRESSION_DIR / "splice_junctions" / "introner_boundary_support"
MICROC_INSULATION_DIR = MICROC_DIR / "qc" / "insulation"
MICROC_TRACK_DIR = MICROC_DIR / "tracks"
MICROC_REF_DIR = MICROC_DIR / "references"

RF_PRESENT_ABSENT_SAMPLES = [RF_TRAIN_SAMPLE, RF_TEST_SAMPLE]


def rf_present_absent_stem(sample):
    """Return filename stem for a sample's present/absent feature matrix."""
    stem = f"{sample}_present_absent_direct"
    if sample == RF_TEST_SAMPLE and RF_EXCLUDED_CONTIGS.get(sample):
        stem += "_exclude28"
    return f"{stem}.mask{RF_MASK}.{RF_WINDOWS_TAG}"


def rf_present_absent_matrix(sample, suffix=""):
    return RF_FEATURES_DIR / f"{rf_present_absent_stem(sample)}{suffix}.tsv"


def rf_present_absent_summary(sample, suffix=""):
    return RF_FEATURES_DIR / f"{rf_present_absent_stem(sample)}{suffix}.summary.tsv"


def rf_enriched_matrix(sample):
    if RF_INCLUDE_EXPRESSION:
        return rf_present_absent_matrix(sample, ".codon_microc_tpm")
    return rf_present_absent_matrix(sample, ".codon_microc")


RF_PREINSERTION_MATRIX = (
    RF_FEATURES_DIR
    / f"rcc1749_preinsertion_sequence_features.mask{RF_MASK}.orientnorm.exclude28.tsv"
)
RF_PREINSERTION_SUMMARY = (
    RF_FEATURES_DIR
    / f"rcc1749_preinsertion_sequence_features.mask{RF_MASK}.orientnorm.exclude28.summary.tsv"
)

LOFO_OUTPUT_DIR = (
    RF_MODELS_DIR
    / f"leave_one_family_out_{RF_TRAIN_SAMPLE.lower()}_present_absent_mask{RF_MASK}_same_gene_negatives"
)
CROSS_SAMPLE_PA_PREFIX = (
    RF_MODELS_DIR
    / f"{RF_TRAIN_SAMPLE.lower()}_train_{RF_TEST_SAMPLE.lower()}_test_present_absent_orientnorm_exclude28_mask{RF_MASK}"
)
PREINSERTION_SCORING_DIR = (
    RF_MODELS_DIR
    / f"{RF_TRAIN_SAMPLE.lower()}_train_{RF_TEST_SAMPLE.lower()}_preinsertion_scoring_mask{RF_MASK}_orientnorm_exclude28"
)
MATCHED_CONTROLS_DIR = (
    RF_MODELS_DIR
    / f"{RF_TRAIN_SAMPLE.lower()}_train_{RF_TEST_SAMPLE.lower()}_preinsertion_matched_controls_mask{RF_MASK}_orientnorm_exclude28"
)
KMER_PERMUTATION_DIR = (
    RF_MODELS_DIR
    / f"{RF_TRAIN_SAMPLE.lower()}_train_{RF_TEST_SAMPLE.lower()}_preinsertion_same_gene_gc_kmer_group_permutation"
)


def rf_exclude_contig_flags(sample):
    contigs = RF_EXCLUDED_CONTIGS.get(sample, [])
    return " ".join(
        f"--exclude-contig {shlex.quote(contig)}" for contig in contigs
    )


def rf_feature_sets_csv(key):
    return ",".join(RF_FEATURE_SETS.get(key, []))


# ============================================================
# FEATURE BUILDING
# ============================================================

rule rf_build_present_absent_ccmp1545:
    """
    Build CCMP1545 present-vs-absent locus feature matrix with boundary support.
    """
    input:
        genotype = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        fa = ASSEMBLIES_DIR / f"{RF_TRAIN_SAMPLE}.vg_paths.fa",
        gtf = lambda wildcards: get_gtf(RF_TRAIN_SAMPLE),
        boundary_support = BOUNDARY_SUPPORT_DIR / f"{RF_TRAIN_SAMPLE}.per_locus.tsv",
    output:
        matrix = rf_present_absent_matrix(RF_TRAIN_SAMPLE),
        summary = rf_present_absent_summary(RF_TRAIN_SAMPLE),
    log:
        RF_LOG_DIR / f"build_present_absent_{RF_TRAIN_SAMPLE}.log",
    params:
        exclude_flags = rf_exclude_contig_flags(RF_TRAIN_SAMPLE),
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {RF_FEATURES_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/build_ccmp1545_present_absent_features.py \
            --genotype-matrix {input.genotype} \
            --fa {input.fa} \
            --gtf {input.gtf} \
            --sample {RF_TRAIN_SAMPLE} \
            --windows {RF_WINDOWS_STR} \
            --junction-mask-bp {RF_MASK} \
            --boundary-max-abs-delta {RF_BOUNDARY_MAX_DELTA} \
            --introner-boundary-support {input.boundary_support} \
            {params.exclude_flags} \
            --output {output.matrix} \
            --summary {output.summary} \
            > {log} 2>&1
        """


rule rf_build_present_absent_rcc1749:
    """
    Build RCC1749 present-vs-absent matrix with orientation normalization.
    """
    input:
        genotype = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        fa = ASSEMBLIES_DIR / f"{RF_TEST_SAMPLE}.vg_paths.fa",
        gtf = lambda wildcards: get_gtf(RF_TEST_SAMPLE),
        mummer_coords = MUMMER_COORDS,
    output:
        matrix = rf_present_absent_matrix(RF_TEST_SAMPLE),
        summary = rf_present_absent_summary(RF_TEST_SAMPLE),
    log:
        RF_LOG_DIR / f"build_present_absent_{RF_TEST_SAMPLE}.log",
    params:
        exclude_flags = rf_exclude_contig_flags(RF_TEST_SAMPLE),
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {RF_FEATURES_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/build_ccmp1545_present_absent_features.py \
            --genotype-matrix {input.genotype} \
            --fa {input.fa} \
            --gtf {input.gtf} \
            --sample {RF_TEST_SAMPLE} \
            --windows {RF_WINDOWS_STR} \
            --junction-mask-bp {RF_MASK} \
            --mummer-coords {input.mummer_coords} \
            --orientation-threshold {RF_ORIENTATION_THRESHOLD} \
            {params.exclude_flags} \
            --output {output.matrix} \
            --summary {output.summary} \
            > {log} 2>&1
        """


rule rf_add_microc_ccmp1545:
    """
    Add Micro-C features to the CCMP1545 present/absent matrix.
    """
    input:
        matrix = rf_present_absent_matrix(RF_TRAIN_SAMPLE),
        insulation = MICROC_INSULATION_DIR / f"{RF_TRAIN_SAMPLE}.{RF_MICROC_PREFIX}.insulation.tsv",
        boundaries = MICROC_INSULATION_DIR / f"{RF_TRAIN_SAMPLE}.{RF_MICROC_PREFIX}.boundaries.bed",
        domains = MICROC_INSULATION_DIR / f"{RF_TRAIN_SAMPLE}.{RF_MICROC_PREFIX}.domains.bed",
        anchors = MICROC_TRACK_DIR / f"{RF_TRAIN_SAMPLE}.anchors.cpm.bedgraph.gz",
        coverage = MICROC_TRACK_DIR / f"{RF_TRAIN_SAMPLE}.coverage.cpm.bedgraph.gz",
        chromsizes = MICROC_REF_DIR / f"{RF_TRAIN_SAMPLE}.chrom.sizes",
    output:
        matrix = rf_present_absent_matrix(RF_TRAIN_SAMPLE, ".microc"),
        summary = RF_FEATURES_DIR / f"{rf_present_absent_stem(RF_TRAIN_SAMPLE)}.microc.join_summary.tsv",
    log:
        RF_LOG_DIR / f"add_microc_{RF_TRAIN_SAMPLE}.log",
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {RF_FEATURES_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/add_microc_features.py \
            --matrix {input.matrix} \
            --microc-dir {MICROC_DIR} \
            --sample {RF_TRAIN_SAMPLE} \
            --resolution {RF_MICROC_RES} \
            --window {RF_MICROC_WINDOW} \
            --flank-windows {RF_MICROC_FLANK_STR} \
            --output {output.matrix} \
            --summary {output.summary} \
            > {log} 2>&1
        """


rule rf_add_microc_rcc1749:
    """
    Add orientation-aware Micro-C features to the RCC1749 present/absent matrix.
    """
    input:
        matrix = rf_present_absent_matrix(RF_TEST_SAMPLE),
        insulation = MICROC_INSULATION_DIR / f"{RF_TEST_SAMPLE}.{RF_MICROC_PREFIX}.insulation.tsv",
        boundaries = MICROC_INSULATION_DIR / f"{RF_TEST_SAMPLE}.{RF_MICROC_PREFIX}.boundaries.bed",
        domains = MICROC_INSULATION_DIR / f"{RF_TEST_SAMPLE}.{RF_MICROC_PREFIX}.domains.bed",
        anchors = MICROC_TRACK_DIR / f"{RF_TEST_SAMPLE}.anchors.cpm.bedgraph.gz",
        coverage = MICROC_TRACK_DIR / f"{RF_TEST_SAMPLE}.coverage.cpm.bedgraph.gz",
        chromsizes = MICROC_REF_DIR / f"{RF_TEST_SAMPLE}.chrom.sizes",
        mummer_coords = MUMMER_COORDS,
    output:
        matrix = rf_present_absent_matrix(RF_TEST_SAMPLE, ".microc"),
        summary = RF_FEATURES_DIR / f"{rf_present_absent_stem(RF_TEST_SAMPLE)}.microc.join_summary.tsv",
    log:
        RF_LOG_DIR / f"add_microc_{RF_TEST_SAMPLE}.log",
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {RF_FEATURES_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/add_microc_features.py \
            --matrix {input.matrix} \
            --microc-dir {MICROC_DIR} \
            --sample {RF_TEST_SAMPLE} \
            --resolution {RF_MICROC_RES} \
            --window {RF_MICROC_WINDOW} \
            --flank-windows {RF_MICROC_FLANK_STR} \
            --mummer-coords {input.mummer_coords} \
            --orientation-threshold {RF_ORIENTATION_THRESHOLD} \
            --output {output.matrix} \
            --summary {output.summary} \
            > {log} 2>&1
        """


rule rf_add_codon_ccmp1545:
    """
    Add codon boundary features to the CCMP1545 Micro-C-enriched matrix.
    """
    input:
        matrix = rf_present_absent_matrix(RF_TRAIN_SAMPLE, ".microc"),
        gtf = lambda wildcards: get_gtf(RF_TRAIN_SAMPLE),
        fa = ASSEMBLIES_DIR / f"{RF_TRAIN_SAMPLE}.vg_paths.fa",
    output:
        matrix = rf_present_absent_matrix(RF_TRAIN_SAMPLE, ".codon_microc"),
        summary = RF_FEATURES_DIR / f"{rf_present_absent_stem(RF_TRAIN_SAMPLE)}.codon_microc.join_summary.tsv",
        unmatched = RF_FEATURES_DIR / f"{rf_present_absent_stem(RF_TRAIN_SAMPLE)}.codon_microc.unmatched.tsv",
    log:
        RF_LOG_DIR / f"add_codon_{RF_TRAIN_SAMPLE}.log",
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {RF_FEATURES_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/add_codon_boundary_features.py \
            --matrix {input.matrix} \
            --gtf {input.gtf} \
            --fa {input.fa} \
            --output {output.matrix} \
            --summary {output.summary} \
            --unmatched {output.unmatched} \
            > {log} 2>&1
        """


rule rf_add_codon_rcc1749:
    """
    Add codon boundary features to the RCC1749 Micro-C-enriched matrix.
    """
    input:
        matrix = rf_present_absent_matrix(RF_TEST_SAMPLE, ".microc"),
        gtf = lambda wildcards: get_gtf(RF_TEST_SAMPLE),
        fa = ASSEMBLIES_DIR / f"{RF_TEST_SAMPLE}.vg_paths.fa",
    output:
        matrix = rf_present_absent_matrix(RF_TEST_SAMPLE, ".codon_microc"),
        summary = RF_FEATURES_DIR / f"{rf_present_absent_stem(RF_TEST_SAMPLE)}.codon_microc.join_summary.tsv",
        unmatched = RF_FEATURES_DIR / f"{rf_present_absent_stem(RF_TEST_SAMPLE)}.codon_microc.unmatched.tsv",
    log:
        RF_LOG_DIR / f"add_codon_{RF_TEST_SAMPLE}.log",
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {RF_FEATURES_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/add_codon_boundary_features.py \
            --matrix {input.matrix} \
            --gtf {input.gtf} \
            --fa {input.fa} \
            --output {output.matrix} \
            --summary {output.summary} \
            --unmatched {output.unmatched} \
            > {log} 2>&1
        """


rule rf_add_expression_ccmp1545:
    """
    Add gene TPM features to the CCMP1545 codon-enriched matrix.
    """
    input:
        matrix = rf_present_absent_matrix(RF_TRAIN_SAMPLE, ".codon_microc"),
        counts = EXPRESSION_DIR / "counts" / "merged_counts_matrix.csv",
        boundary_support = BOUNDARY_SUPPORT_DIR / f"{RF_TRAIN_SAMPLE}.per_locus.tsv",
    output:
        matrix = rf_present_absent_matrix(RF_TRAIN_SAMPLE, ".codon_microc_tpm"),
        summary = RF_FEATURES_DIR / f"{rf_present_absent_stem(RF_TRAIN_SAMPLE)}.codon_microc_tpm.join_summary.tsv",
    log:
        RF_LOG_DIR / f"add_expression_{RF_TRAIN_SAMPLE}.log",
    params:
        skip_psi_flag = "" if RF_INCLUDE_PSI else "--skip-psi",
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {RF_FEATURES_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/add_expression_features.py \
            --matrix {input.matrix} \
            --counts {input.counts} \
            --boundary-support {input.boundary_support} \
            --sample-prefix {RF_EXPR_PREFIX} \
            {params.skip_psi_flag} \
            --output {output.matrix} \
            --summary {output.summary} \
            > {log} 2>&1
        """


rule rf_add_expression_rcc1749:
    """
    Add gene TPM features to the RCC1749 codon-enriched matrix.
    """
    input:
        matrix = rf_present_absent_matrix(RF_TEST_SAMPLE, ".codon_microc"),
        counts = EXPRESSION_DIR / "counts" / "merged_counts_matrix.csv",
        boundary_support = BOUNDARY_SUPPORT_DIR / f"{RF_TRAIN_SAMPLE}.per_locus.tsv",
    output:
        matrix = rf_present_absent_matrix(RF_TEST_SAMPLE, ".codon_microc_tpm"),
        summary = RF_FEATURES_DIR / f"{rf_present_absent_stem(RF_TEST_SAMPLE)}.codon_microc_tpm.join_summary.tsv",
    log:
        RF_LOG_DIR / f"add_expression_{RF_TEST_SAMPLE}.log",
    params:
        skip_psi_flag = "" if RF_INCLUDE_PSI else "--skip-psi",
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {RF_FEATURES_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/add_expression_features.py \
            --matrix {input.matrix} \
            --counts {input.counts} \
            --boundary-support {input.boundary_support} \
            --sample-prefix {RF_EXPR_PREFIX} \
            {params.skip_psi_flag} \
            --output {output.matrix} \
            --summary {output.summary} \
            > {log} 2>&1
        """


rule rf_build_preinsertion_features:
    """
    Build RCC1749 empty-ortholog pre-insertion sequence feature matrix.
    """
    input:
        genotype = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        gtf = lambda wildcards: get_gtf(RF_TEST_SAMPLE),
        fa = ASSEMBLIES_DIR / f"{RF_TEST_SAMPLE}.vg_paths.fa",
        mummer_coords = MUMMER_COORDS,
    output:
        matrix = RF_PREINSERTION_MATRIX,
        summary = RF_PREINSERTION_SUMMARY,
    log:
        RF_LOG_DIR / "build_preinsertion_features.log",
    params:
        exclude_flags = rf_exclude_contig_flags(RF_TEST_SAMPLE),
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {RF_FEATURES_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/build_preinsertion_features.py \
            --genotype-matrix {input.genotype} \
            --gtf {input.gtf} \
            --fa {input.fa} \
            --sample {RF_TEST_SAMPLE} \
            --reference-sample {RF_TRAIN_SAMPLE} \
            --windows {RF_WINDOWS_STR} \
            --junction-mask-bp {RF_MASK} \
            --controls-per-positive {RF_MODELING[preinsertion_controls_per_positive]} \
            --mummer-coords {input.mummer_coords} \
            --orientation-threshold {RF_ORIENTATION_THRESHOLD} \
            {params.exclude_flags} \
            --output {output.matrix} \
            --summary {output.summary} \
            > {log} 2>&1
        """


# ============================================================
# MODELING
# ============================================================

rule rf_leave_one_family_out:
    """
    Leave-one-introner-family-out CV within the training sample.
    """
    input:
        matrix = rf_enriched_matrix(RF_TRAIN_SAMPLE),
    output:
        summary = LOFO_OUTPUT_DIR / "summary.tsv",
        predictions = LOFO_OUTPUT_DIR / "predictions.tsv",
        features = LOFO_OUTPUT_DIR / "features.tsv",
        skipped = LOFO_OUTPUT_DIR / "skipped_families.tsv",
    log:
        RF_LOG_DIR / "leave_one_family_out.log",
    params:
        feature_sets = rf_feature_sets_csv("within_sample"),
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {LOFO_OUTPUT_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/run_leave_one_family_out.py \
            --matrix {input.matrix} \
            --windows {RF_WINDOWS_STR} \
            --feature-sets {params.feature_sets} \
            --negative-mode {RF_MODELING[lofo_negative_mode]} \
            --min-test-pos {RF_MODELING[lofo_min_test_pos]} \
            --min-test-neg {RF_MODELING[lofo_min_test_neg]} \
            --n-estimators {RF_MODELING[n_estimators]} \
            --random-state {RF_MODELING[random_state]} \
            --output-dir {LOFO_OUTPUT_DIR} \
            > {log} 2>&1
        """


rule rf_cross_sample_present_absent:
    """
    Train on CCMP1545 present/absent; evaluate on RCC1749 present/absent.
    """
    input:
        train_matrix = rf_enriched_matrix(RF_TRAIN_SAMPLE),
        test_matrix = rf_enriched_matrix(RF_TEST_SAMPLE),
    output:
        summary = f"{CROSS_SAMPLE_PA_PREFIX}.summary.tsv",
        predictions = f"{CROSS_SAMPLE_PA_PREFIX}.predictions.tsv",
        features = f"{CROSS_SAMPLE_PA_PREFIX}.features.tsv",
    log:
        RF_LOG_DIR / "cross_sample_present_absent.log",
    params:
        feature_sets = rf_feature_sets_csv("cross_sample_transfer"),
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {RF_MODELS_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/run_cross_sample_present_absent.py \
            --train-matrix {input.train_matrix} \
            --test-matrix {input.test_matrix} \
            --windows {RF_WINDOWS_STR} \
            --feature-sets {params.feature_sets} \
            --n-estimators {RF_MODELING[n_estimators]} \
            --random-state {RF_MODELING[random_state]} \
            --output-prefix {CROSS_SAMPLE_PA_PREFIX} \
            > {log} 2>&1
        """


rule rf_score_preinsertion:
    """
    Train CCMP1545 present/absent RF and score RCC1749 pre-insertion sites.
    """
    input:
        train_matrix = rf_enriched_matrix(RF_TRAIN_SAMPLE),
        test_matrix = RF_PREINSERTION_MATRIX,
    output:
        summary = PREINSERTION_SCORING_DIR / "summary.tsv",
        predictions = PREINSERTION_SCORING_DIR / "predictions.tsv",
        features = PREINSERTION_SCORING_DIR / "features.tsv",
        config = PREINSERTION_SCORING_DIR / "config.txt",
    log:
        RF_LOG_DIR / "score_preinsertion.log",
    params:
        feature_sets = rf_feature_sets_csv("cross_sample_sequence"),
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {PREINSERTION_SCORING_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/run_ccmp_model_score_preinsertion.py \
            --train-matrix {input.train_matrix} \
            --test-matrix {input.test_matrix} \
            --windows {RF_WINDOWS_STR} \
            --feature-sets {params.feature_sets} \
            --n-estimators {RF_MODELING[n_estimators]} \
            --random-state {RF_MODELING[random_state]} \
            --output-dir {PREINSERTION_SCORING_DIR} \
            > {log} 2>&1
        """


rule rf_preinsertion_matched_control_scoring:
    """
    Score empty orthologs against GC-matched exonic controls.
    """
    input:
        positive_matrix = RF_PREINSERTION_MATRIX,
        train_matrix = rf_enriched_matrix(RF_TRAIN_SAMPLE),
        genotype = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        gtf = lambda wildcards: get_gtf(RF_TEST_SAMPLE),
        fa = ASSEMBLIES_DIR / f"{RF_TEST_SAMPLE}.vg_paths.fa",
        mummer_coords = MUMMER_COORDS,
    output:
        summary = MATCHED_CONTROLS_DIR / "summary.tsv",
        control_summary = MATCHED_CONTROLS_DIR / "control_candidate_summary.tsv",
        gc_predictions = MATCHED_CONTROLS_DIR / "gc_only.predictions.tsv",
        same_gene_predictions = MATCHED_CONTROLS_DIR / "same_gene_gc.predictions.tsv",
        same_gene_matched = MATCHED_CONTROLS_DIR / "same_gene_gc.matched_dataset.tsv",
    log:
        RF_LOG_DIR / "preinsertion_matched_controls.log",
    params:
        feature_sets = rf_feature_sets_csv("cross_sample_sequence"),
        exclude_flags = rf_exclude_contig_flags(RF_TEST_SAMPLE),
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {MATCHED_CONTROLS_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/run_preinsertion_matched_control_scoring.py \
            --positive-matrix {input.positive_matrix} \
            --train-matrix {input.train_matrix} \
            --genotype-matrix {input.genotype} \
            --gtf {input.gtf} \
            --fa {input.fa} \
            --mummer-coords {input.mummer_coords} \
            --windows {RF_WINDOWS_STR} \
            --gc-col {RF_MODELING[matched_gc_col]} \
            --gc-tolerance {RF_MODELING[matched_gc_tolerance]} \
            --feature-sets {params.feature_sets} \
            --n-estimators {RF_MODELING[n_estimators]} \
            --random-state {RF_MODELING[random_state]} \
            {params.exclude_flags} \
            --output-dir {MATCHED_CONTROLS_DIR} \
            > {log} 2>&1
        """


rule rf_cross_sample_kmer_group_permutation:
    """
    Grouped k-mer permutation importance on cross-sample matched test set.
    """
    input:
        train_matrix = rf_enriched_matrix(RF_TRAIN_SAMPLE),
        test_matrix = MATCHED_CONTROLS_DIR / "same_gene_gc.matched_dataset.tsv",
    output:
        baseline = KMER_PERMUTATION_DIR / "baseline.tsv",
        importance = KMER_PERMUTATION_DIR / "kmer_permutation_importance.tsv",
        importance_by_repeat = KMER_PERMUTATION_DIR / "kmer_permutation_importance_by_repeat.tsv",
        feature_groups = KMER_PERMUTATION_DIR / "feature_groups.tsv",
        baseline_predictions = KMER_PERMUTATION_DIR / "baseline_predictions.tsv",
    log:
        RF_LOG_DIR / "cross_sample_kmer_permutation.log",
    conda: "../envs/random_forest.yaml"
    shell:
        """
        mkdir -p {KMER_PERMUTATION_DIR} {RF_LOG_DIR}
        export PYTHONPATH={RF_SCRIPTS}:${{PYTHONPATH:-}}

        python {RF_SCRIPTS}/run_cross_sample_kmer_group_permutation.py \
            --train-matrix {input.train_matrix} \
            --test-matrix {input.test_matrix} \
            --feature-set kmer_left_right \
            --grouping kmer \
            --n-repeats {RF_MODELING[kmer_permutation_repeats]} \
            --n-estimators {RF_MODELING[n_estimators]} \
            --random-state {RF_MODELING[random_state]} \
            --output-dir {KMER_PERMUTATION_DIR} \
            > {log} 2>&1
        """


# ============================================================
# TARGET RULE
# ============================================================

rule random_forest_all:
    """
    Run the full random forest classifier analysis pipeline.
    """
    input:
        rf_enriched_matrix(RF_TRAIN_SAMPLE),
        rf_enriched_matrix(RF_TEST_SAMPLE),
        RF_PREINSERTION_MATRIX,
        LOFO_OUTPUT_DIR / "summary.tsv",
        f"{CROSS_SAMPLE_PA_PREFIX}.summary.tsv",
        PREINSERTION_SCORING_DIR / "summary.tsv",
        MATCHED_CONTROLS_DIR / "summary.tsv",
        KMER_PERMUTATION_DIR / "kmer_permutation_importance.tsv",
