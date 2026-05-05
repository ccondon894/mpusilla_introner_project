# ============================================================================
# 16_sfs_inference.smk
# Tree-free SFS-based inference of introner gain/loss dynamics
#
# Fits a single-population demographic model to the 4D-site SFS, then jointly
# estimates per-locus-class gain rate (theta_lambda) and loss rate (theta_mu)
# from the introner two-state SFS as a Poisson-Random-Field gain+loss mixture.
# Compares introners against non-introner introns (the "control" locus class)
# under matched samples and matched demographic background.
#
# Pipeline structure:
#
#   per subset (e.g. full_g1, cluster_0, cluster_1, inner_quartet):
#     1. tabulate folded 4D SFS from upstream all-sites 4D VCF
#     2. fit demography (constant + 2/3-epoch + growth + bottlegrowth) by AIC
#     3. validate demography (Tajima's D / pi / theta_W vs Poisson sims)
#
#     for each locus class (introner, non_introner_intron):
#       4. tabulate locus two-state SFS from genotype matrix
#       5. fit Poisson gain/loss MLE using fitted demography
#       6. (optional) bin-level monotonicity-constrained fit
#
#   aggregate:
#     7. compare introners vs non-introner introns side-by-side
#
# Optional one-off diagnostics (not in default DAG):
#     - check_g1_substructure  (PCA + F_ST on G1 4D SNPs)
#     - pairwise_ibs           (within-cluster pairwise IBS)
#
# Dependencies:
#   - results/genotyping/genotype_matrix.final.tsv (introner matrix)
#   - results/evolution/non_introner_introns/intron_genotype_matrix.tsv
#       (built by 13_non_introner_introns.smk)
#   - 4D all-sites VCF for Group 1 (path configurable; default is the upstream
#     introner-genotyping-pipeline output — see config["sfs_inference"]["fourd_vcf"])
#
# Conda env: workflow/envs/moments.yaml
# ============================================================================

SFS_DIR = EVOLUTION_DIR / "sfs_inference"
SFS_SCRIPTS = str(PROJECT_ROOT / "scripts" / "evolution" / "sfs_inference")

# 4D VCF: defaults to the upstream introner-genotyping-pipeline output unless
# overridden in config. Path is plain gzip (not bgzip), 11 G1 samples incl.
# CCMP1545, all-sites GVCF restricted to 4-fold-degenerate sites with MT
# excluded. See `RESULTS_SUMMARY.md` for upstream provenance.
FOURD_VCF = config.get("sfs_inference", {}).get(
    "fourd_vcf",
    "/scratch1/chris/introner-genotyping-pipeline/popgen_snp_data/"
    "mpusilla.snps.4d.notMT.group1.allsites.vcf.gz",
)

# Sample subsets: full G1 plus the substructure-corrected sub-clusters
# discovered by PCA + IBS (see check_g1_substructure / pairwise_ibs rules).
# Override in config under sfs_inference.subsets if needed.
DEFAULT_SFS_SUBSETS = {
    "full_g1": GROUP1_SAMPLES,
    "cluster_0": ["RCC114", "RCC2482", "RCC629", "RCC693", "RCC833"],
    "cluster_1": ["CCMP1545", "RCC1614", "RCC1698", "RCC373", "RCC465", "RCC692"],
    "inner_quartet": ["RCC1614", "RCC1698", "RCC373", "RCC465"],
}
SFS_SUBSETS = config.get("sfs_inference", {}).get("subsets", DEFAULT_SFS_SUBSETS)

LOCUS_CLASSES = ["introner", "non_introner_intron"]

# Bin-level model bootstrap settings
BIN_LEVEL_BOOTSTRAP_REPS = config.get("sfs_inference", {}).get("bootstrap_reps", 1000)
BIN_LEVEL_SEED = config.get("sfs_inference", {}).get("seed", 42)

# Wildcard constraints (subset names + locus class names)
wildcard_constraints:
    sfs_subset = "|".join(SFS_SUBSETS.keys()),
    locus_class = "|".join(LOCUS_CLASSES),


# ============================================================================
# 4D SFS + demography (per subset)
# ============================================================================

rule sfs_tabulate_4d:
    """Tabulate folded 4D SFS for a sample subset from the upstream all-sites VCF."""
    input:
        vcf = FOURD_VCF,
    output:
        npy = SFS_DIR / "{sfs_subset}" / "4d_sfs.npy",
        plot = SFS_DIR / "{sfs_subset}" / "4d_sfs.png",
        summary = SFS_DIR / "{sfs_subset}" / "4d_sfs_summary.txt",
    params:
        samples = lambda wc: ",".join(SFS_SUBSETS[wc.sfs_subset]),
    conda:
        "../envs/moments.yaml"
    shell:
        """
        mkdir -p $(dirname {output.npy}) && \\
        python {SFS_SCRIPTS}/tabulate_4d_sfs.py \\
            --vcf {input.vcf} \\
            --samples {params.samples} \\
            --output {output.npy} \\
            --plot {output.plot} \\
            --summary {output.summary}
        """


rule sfs_fit_demography:
    """Fit constant + two/three-epoch + growth + bottlegrowth models, pick best by AIC.

    Outputs a per-bin Pearson residuals TSV alongside the AIC table — these
    residuals are how we originally diagnosed within-G1 substructure.
    """
    input:
        sfs = SFS_DIR / "{sfs_subset}" / "4d_sfs.npy",
    output:
        aic = SFS_DIR / "{sfs_subset}" / "demography_aic.tsv",
        best = SFS_DIR / "{sfs_subset}" / "demography_best.pkl",
        plot = SFS_DIR / "{sfs_subset}" / "demography_diagnostic.png",
        residuals = SFS_DIR / "{sfs_subset}" / "demography_residuals.tsv",
    conda:
        "../envs/moments.yaml"
    shell:
        """
        python {SFS_SCRIPTS}/fit_demography.py \\
            --sfs {input.sfs} \\
            --aic {output.aic} \\
            --best {output.best} \\
            --plot {output.plot} \\
            --residuals {output.residuals}
        """


rule sfs_validate_demography:
    """Compare observed Tajima's D / pi / theta_W to Poisson-resampled simulations
    under the fitted demography. Underestimates true coalescent variance, so a
    rejection here is conservative."""
    input:
        sfs = SFS_DIR / "{sfs_subset}" / "4d_sfs.npy",
        demog = SFS_DIR / "{sfs_subset}" / "demography_best.pkl",
    output:
        summary = SFS_DIR / "{sfs_subset}" / "demography_validation.tsv",
        plot = SFS_DIR / "{sfs_subset}" / "demography_validation.png",
    conda:
        "../envs/moments.yaml"
    shell:
        """
        python {SFS_SCRIPTS}/validate_demography.py \\
            --sfs {input.sfs} \\
            --demog {input.demog} \\
            --summary {output.summary} \\
            --plot {output.plot}
        """


# ============================================================================
# Locus-class SFS (per subset, per locus class)
# ============================================================================

rule sfs_tabulate_introner:
    """Tabulate the introner two-state SFS from the genotype matrix.

    Applies the within_group_status ∈ {consistent, singleton} filter to
    exclude low_identity orthologs (where parallel insertion is plausible)
    and drops loci with any presence=3 (missing) call in the subset.
    """
    input:
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
    output:
        sfs = SFS_DIR / "{sfs_subset}" / "introner" / "sfs.npy",
        summary = SFS_DIR / "{sfs_subset}" / "introner" / "sfs_summary.tsv",
        plot = SFS_DIR / "{sfs_subset}" / "introner" / "sfs.png",
    params:
        samples = lambda wc: ",".join(SFS_SUBSETS[wc.sfs_subset]),
    conda:
        "../envs/moments.yaml"
    shell:
        """
        mkdir -p $(dirname {output.sfs}) && \\
        python {SFS_SCRIPTS}/tabulate_introner_sfs.py \\
            --matrix {input.matrix} \\
            --samples {params.samples} \\
            --output {output.sfs} \\
            --summary {output.summary} \\
            --plot {output.plot}
        """


rule sfs_tabulate_non_introner_intron:
    """Tabulate the non-introner intron two-state SFS — the locus-class control."""
    input:
        matrix = EVOLUTION_DIR / "non_introner_introns" / "intron_genotype_matrix.tsv",
    output:
        sfs = SFS_DIR / "{sfs_subset}" / "non_introner_intron" / "sfs.npy",
        summary = SFS_DIR / "{sfs_subset}" / "non_introner_intron" / "sfs_summary.tsv",
        plot = SFS_DIR / "{sfs_subset}" / "non_introner_intron" / "sfs.png",
    params:
        samples = lambda wc: ",".join(SFS_SUBSETS[wc.sfs_subset]),
    conda:
        "../envs/moments.yaml"
    shell:
        """
        mkdir -p $(dirname {output.sfs}) && \\
        python {SFS_SCRIPTS}/tabulate_non_introner_sfs.py \\
            --matrix {input.matrix} \\
            --samples {params.samples} \\
            --output {output.sfs} \\
            --summary {output.summary} \\
            --plot {output.plot}
        """


# ============================================================================
# Gain/loss MLE (per subset, per locus class)
# ============================================================================

rule sfs_fit_gain_loss:
    """Joint MLE of (theta_lambda, theta_mu) under the fitted demographic SFS shape.

    Output: per-bin observed/expected, G-test goodness-of-fit, and a fit plot.
    """
    input:
        sfs = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "sfs.npy",
        summary = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "sfs_summary.tsv",
        demog = SFS_DIR / "{sfs_subset}" / "demography_best.pkl",
    output:
        estimates = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "theta_estimates.tsv",
        gof = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "gain_loss_gof.tsv",
        plot = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "gain_loss_fit.png",
    conda:
        "../envs/moments.yaml"
    shell:
        """
        python {SFS_SCRIPTS}/fit_gain_loss.py \\
            --sfs {input.sfs} \\
            --summary {input.summary} \\
            --demog {input.demog} \\
            --estimates {output.estimates} \\
            --gof {output.gof} \\
            --plot {output.plot}
        """


rule sfs_fit_bin_level:
    """Bin-level (monotonicity-constrained) gain/loss MLE.

    Alternative to the parametric model: relaxes the demography-conditional
    SFS shape to per-bin alpha (gain) and beta (loss) under monotonicity
    constraints. Reports an LP-based identification interval for the ratio
    when the saturated MLE is reachable (the model is non-identified there).
    """
    input:
        sfs = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "sfs.npy",
        summary = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "sfs_summary.tsv",
        parametric = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "theta_estimates.tsv",
    output:
        pkl = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "bin_level" / "fit.pkl",
        tsv = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "bin_level" / "fit.tsv",
        json = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "bin_level" / "summary.json",
        plot = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "bin_level" / "residuals.png",
        warnings = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "bin_level" / "warnings.txt",
        comparison = SFS_DIR / "{sfs_subset}" / "{locus_class}" / "bin_level" / "comparison_with_parametric.tsv",
    params:
        bootstrap_reps = BIN_LEVEL_BOOTSTRAP_REPS,
        seed = BIN_LEVEL_SEED,
    conda:
        "../envs/moments.yaml"
    shell:
        """
        mkdir -p $(dirname {output.pkl}) && \\
        python {SFS_SCRIPTS}/fit_bin_level_model.py \\
            --sfs {input.sfs} \\
            --summary {input.summary} \\
            --parametric_estimates {input.parametric} \\
            --out_pkl {output.pkl} \\
            --out_tsv {output.tsv} \\
            --out_summary {output.json} \\
            --out_plot {output.plot} \\
            --out_warnings {output.warnings} \\
            --out_comparison {output.comparison} \\
            --bootstrap_reps {params.bootstrap_reps} \\
            --seed {params.seed}
        """


# ============================================================================
# Cross-locus-class comparison
# ============================================================================

rule sfs_compare_introner_vs_non_introner:
    """Side-by-side introner vs non-introner intron comparison across subsets.

    Aggregates all per-(subset, locus_class) fits and produces:
      - results/evolution/sfs_inference/comparison/introner_vs_non_introner.tsv
      - results/figures/introner_vs_non_introner_sfs.png
    """
    input:
        introner_estimates = expand(
            str(SFS_DIR / "{sfs_subset}" / "introner" / "theta_estimates.tsv"),
            sfs_subset=list(SFS_SUBSETS.keys()),
        ),
        non_introner_estimates = expand(
            str(SFS_DIR / "{sfs_subset}" / "non_introner_intron" / "theta_estimates.tsv"),
            sfs_subset=list(SFS_SUBSETS.keys()),
        ),
        introner_sfs = expand(
            str(SFS_DIR / "{sfs_subset}" / "introner" / "sfs.npy"),
            sfs_subset=list(SFS_SUBSETS.keys()),
        ),
        non_introner_sfs = expand(
            str(SFS_DIR / "{sfs_subset}" / "non_introner_intron" / "sfs.npy"),
            sfs_subset=list(SFS_SUBSETS.keys()),
        ),
        introner_summaries = expand(
            str(SFS_DIR / "{sfs_subset}" / "introner" / "sfs_summary.tsv"),
            sfs_subset=list(SFS_SUBSETS.keys()),
        ),
        non_introner_summaries = expand(
            str(SFS_DIR / "{sfs_subset}" / "non_introner_intron" / "sfs_summary.tsv"),
            sfs_subset=list(SFS_SUBSETS.keys()),
        ),
        introner_gofs = expand(
            str(SFS_DIR / "{sfs_subset}" / "introner" / "gain_loss_gof.tsv"),
            sfs_subset=list(SFS_SUBSETS.keys()),
        ),
        non_introner_gofs = expand(
            str(SFS_DIR / "{sfs_subset}" / "non_introner_intron" / "gain_loss_gof.tsv"),
            sfs_subset=list(SFS_SUBSETS.keys()),
        ),
    output:
        tsv = SFS_DIR / "comparison" / "introner_vs_non_introner.tsv",
        png = FIGURES_DIR / "introner_vs_non_introner_sfs.png",
    params:
        root = lambda wc, output: str(Path(output.tsv).parent.parent),
        subsets = ",".join(SFS_SUBSETS.keys()),
    conda:
        "../envs/moments.yaml"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) {FIGURES_DIR} && \\
        python {SFS_SCRIPTS}/compare_introner_vs_non_introner.py \\
            --root {params.root} \\
            --subsets {params.subsets} \\
            --out_tsv {output.tsv} \\
            --out_png {output.png}
        """


# ============================================================================
# Optional / one-off diagnostic rules — not in default DAG
# ============================================================================

rule sfs_check_g1_substructure:
    """One-off diagnostic: PCA + Hudson F_ST on G1 4D segregating SNPs.

    Used to discover the G1 substructure (5/6 cluster_0 vs cluster_1 split,
    F_ST = 0.62) that motivates the cluster-based SFS subsets above.
    Run manually:
        snakemake --use-conda results/evolution/sfs_inference/g1_substructure/g1_pca_fst.png
    """
    input:
        vcf = FOURD_VCF,
    output:
        pca_plot = SFS_DIR / "g1_substructure" / "g1_substructure.png",
        pca_table = SFS_DIR / "g1_substructure" / "g1_pca.tsv",
        fst = SFS_DIR / "g1_substructure" / "g1_fst.tsv",
        sfs_decomp_plot = SFS_DIR / "g1_substructure" / "g1_sfs_decomp.png",
        sfs_decomp_table = SFS_DIR / "g1_substructure" / "g1_sfs_decomp.tsv",
        summary_plot = SFS_DIR / "g1_substructure" / "g1_pca_fst.png",
    params:
        samples = ",".join(GROUP1_SAMPLES),
    conda:
        "../envs/moments.yaml"
    shell:
        """
        mkdir -p $(dirname {output.pca_plot}) && \\
        python {SFS_SCRIPTS}/check_g1_substructure.py \\
            --vcf {input.vcf} \\
            --samples {params.samples} \\
            --pca_plot {output.pca_plot} \\
            --pca_table {output.pca_table} \\
            --fst {output.fst} \\
            --sfs_decomp_plot {output.sfs_decomp_plot} \\
            --sfs_decomp_table {output.sfs_decomp_table} \\
            --summary_plot {output.summary_plot}
        """


rule sfs_pairwise_ibs:
    """One-off diagnostic: pairwise IBS within a subset (e.g. cluster_1 to find
    the inner-quartet vs outgroup-pair split). Run manually:
        snakemake --use-conda results/evolution/sfs_inference/cluster_1/ibs_diagnostic.png
    """
    input:
        vcf = FOURD_VCF,
    output:
        matrix = SFS_DIR / "{sfs_subset}" / "ibs_matrix.tsv",
        pairs = SFS_DIR / "{sfs_subset}" / "ibs_pairs.tsv",
        bin2_decomp = SFS_DIR / "{sfs_subset}" / "ibs_bin2_ownership.tsv",
        plot = SFS_DIR / "{sfs_subset}" / "ibs_diagnostic.png",
    params:
        samples = lambda wc: ",".join(SFS_SUBSETS[wc.sfs_subset]),
    conda:
        "../envs/moments.yaml"
    shell:
        """
        python {SFS_SCRIPTS}/pairwise_ibs.py \\
            --vcf {input.vcf} \\
            --samples {params.samples} \\
            --ibs_matrix {output.matrix} \\
            --ibs_summary {output.pairs} \\
            --bin2_decomp {output.bin2_decomp} \\
            --plot {output.plot}
        """


# ============================================================================
# Master target
# ============================================================================

rule sfs_inference_all:
    """Build the full SFS inference pipeline: per-subset demography + per-class
    SFS + gain/loss fits + comparison."""
    input:
        comparison_tsv = SFS_DIR / "comparison" / "introner_vs_non_introner.tsv",
        comparison_png = FIGURES_DIR / "introner_vs_non_introner_sfs.png",
        validations = expand(
            str(SFS_DIR / "{sfs_subset}" / "demography_validation.tsv"),
            sfs_subset=list(SFS_SUBSETS.keys()),
        ),
        bin_level = expand(
            str(SFS_DIR / "{sfs_subset}" / "{locus_class}" / "bin_level" / "summary.json"),
            sfs_subset=list(SFS_SUBSETS.keys()),
            locus_class=LOCUS_CLASSES,
        ),
