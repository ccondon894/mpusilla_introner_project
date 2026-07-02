# ============================================================
# 24_selection.smk - Selection & Population Structure Analysis
# ============================================================
#
# Performs various population genetics analyses:
# - Tajima's D (neutrality test)
# - Site frequency spectra
# - Principal component analysis (PCA)
# - Linkage disequilibrium (LD)
# - Recombination rate comparisons
#
# Adapted from:
# - /scratch1/chris/introner_vis/basic_popgen/
# - /scratch1/chris/introner_vis/pca/
# - /scratch1/chris/introner_vis/ld/
# - /scratch1/chris/introner_vis/polymorphism_analysis/
# - /scratch1/chris/introner_vis/recombination_analysis/
#
# ============================================================

import os
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Output directories
SNP_DIR = RESULTS / "snp_popgen"
VCF_DIR = SNP_DIR / "vcf"
SELECTION_DIR = SNP_DIR / "selection"
PCA_DIR = SELECTION_DIR / "pca"
LD_DIR = SELECTION_DIR / "ld"
SFS_DIR = SELECTION_DIR / "sfs"
RECOMB_DIR = SELECTION_DIR / "recombination"
SWEEP_DIR = SELECTION_DIR / "sweeps"
FIXED_SNP_SWEEP_DIR = SWEEP_DIR / "fixed_snp_control"
SELECTION_LOG_DIR = SNP_DIR / "logs" / "selection"
SNPEFF_DIR = SNP_DIR / "snpeff"

# Pyrho recombination map paths
PYRHO_CCMP1545 = config["paths"]["pyrho"]["ccmp1545"]
PYRHO_RCC1749 = config["paths"]["pyrho"]["rcc1749"]

# Group 1 polymorphic introner sweep screen parameters
SWEEP_CONFIG = config["params"].get("selection_sweeps", {})
SWEEP_WINDOW_SIZES = SWEEP_CONFIG.get("window_sizes", [10000, 25000, 50000])
SWEEP_MIN_CALLABLE_SITES = SWEEP_CONFIG.get("min_callable_sites", 20)
SWEEP_BACKGROUNDS_PER_FOCAL = SWEEP_CONFIG.get("background_windows_per_focal", 50)
SWEEP_FIXED_SNP_MIN_SPACING = SWEEP_CONFIG.get("fixed_snp_min_spacing", max(SWEEP_WINDOW_SIZES))
SWEEP_ACCEPTED_WITHIN_STATUSES = SWEEP_CONFIG.get(
    "accepted_within_statuses", ["consistent", "singleton"]
)

# SweepFinder2 genome-wide CLR scan parameters
SF2_DIR = SWEEP_DIR / "sweepfinder"
SF2_CONFIG = config["params"].get("sweepfinder", {})
SF2_GRID_SPACING = SF2_CONFIG.get("grid_spacing", 1000)
SF2_CLR_PERCENTILE = SF2_CONFIG.get("clr_percentile", 99)
SF2_REGION_MERGE_GAP = SF2_CONFIG.get("region_merge_gap", 2000)
SF2_PERMUTATIONS = SF2_CONFIG.get("permutations", 1000)
SF2_MIN_CALLED_HAPLOTYPES = SF2_CONFIG.get("min_called_haplotypes", 8)
SF2_FLANK_SIZE = SF2_CONFIG.get("flank_size", 100)
SF2_SEED = SF2_CONFIG.get("seed", 42)
CCMP1545_GRAPH_FAI = PROJECT_ROOT / f"{config['paths']['references']['ccmp1545_graph']}.fai"


# ============================================================
# BASIC POPULATION GENETICS STATISTICS
# ============================================================

rule calculate_tajimas_d:
    """
    Calculate Tajima's D statistic for neutrality testing.

    Tajima's D compares two estimates of genetic diversity:
    - Positive D: Balancing selection or population contraction
    - Negative D: Purifying selection or population expansion
    - D ~ 0: Neutral evolution
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.group1.vcf.gz"
    output:
        tsv = SELECTION_DIR / "tajimas_d.tsv",
        summary = SELECTION_DIR / "tajimas_d_summary.txt"
    log:
        SELECTION_LOG_DIR / "tajimas_d.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {SELECTION_DIR}
        mkdir -p {SELECTION_LOG_DIR}

        python {PROJECT_ROOT}/scripts/popgen/basic_popgen/tajima_d.py \
            --vcf {input.vcf} \
            --output {output.tsv} \
            --summary {output.summary} \
            2> {log}
        """


rule plot_sfs_by_class:
    """
    Plots the SFS frequencies for synonymous/nonsynonymous and introner allele frequencies
    """
    input:
        snpeff_vcf = SNPEFF_DIR / "mpusilla.snps.snpEff.no_MT.vcf",
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        non_introner_matrix = RESULTS / "evolution" / "non_introner_introns" / "intron_genotype_matrix.tsv",
        color_guide = PROJECT_ROOT / "master_figure_color_guide.tsv",
        script = PROJECT_ROOT / "scripts" / "popgen" / "polymorphism_analysis" / "plot_afs_by_class.py"
    output:
        tsv = SFS_DIR / "afs_by_class.tsv",
        pdf = FIGURES_DIR / "snp_popgen" / "afs_by_class.pdf",
        png = FIGURES_DIR / "snp_popgen" / "afs_by_class.png"
    params:
        group1_str = ",".join(GROUP1_SAMPLES),
        outgroup_str = ",".join(GROUP2_SAMPLES)
    log:
        SELECTION_LOG_DIR / "afs_by_class.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        python {input.script} \
            --snpeff_vcf {input.snpeff_vcf} \
            --genotype_matrix {input.genotype_matrix} \
            --non_introner_matrix {input.non_introner_matrix} \
            --group1 {params.group1_str} \
            --outgroup {params.outgroup_str} \
            --output_tsv {output.tsv} \
            --output_pdf {output.pdf} \
            --output_png {output.png} \
            --color_guide {input.color_guide} \
            2> {log}
        """


# ============================================================
# GROUP 1 POLYMORPHIC INTRONER SWEEP SCREEN
# ============================================================

rule build_group1_introner_sweep_targets:
    """
    Build CCMP1545-anchored Group 1-only polymorphic introner targets.

    Retains complete Group 1 polymorphisms with accepted within-group orthology,
    no Group 2 presence calls, and no overlap with the CCMP1545 mating-type
    region.
    """
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        tsv = SWEEP_DIR / "group1_polymorphic_introner_targets.tsv",
        bed = SWEEP_DIR / "group1_polymorphic_introner_targets.bed"
    params:
        group1 = ",".join(GROUP1_SAMPLES),
        group2 = ",".join(GROUP2_SAMPLES),
        reference_sample = REFERENCE,
        accepted_statuses = ",".join(SWEEP_ACCEPTED_WITHIN_STATUSES),
        mating_contig = lambda wildcards: f"{REFERENCE}#0#{config['mating_type_region']['scaffold']}",
        mating_start = config["mating_type_region"]["start"],
        mating_end = config["mating_type_region"]["end"]
    log:
        SELECTION_LOG_DIR / "build_group1_introner_sweep_targets.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {SWEEP_DIR} {SELECTION_LOG_DIR}

        python {PROJECT_ROOT}/scripts/popgen/selection_analysis/build_group1_introner_sweep_targets.py \
            --genotype-matrix {input.genotype_matrix} \
            --output-tsv {output.tsv} \
            --output-bed {output.bed} \
            --group1-samples {params.group1} \
            --group2-samples {params.group2} \
            --reference-sample {params.reference_sample} \
            --accepted-within-statuses {params.accepted_statuses} \
            --mating-contig {params.mating_contig} \
            --mating-start {params.mating_start} \
            --mating-end {params.mating_end} \
            > {log} 2>&1
        """


rule calculate_group1_introner_sweep_windows:
    """
    Compare sweep statistics around Group 1 polymorphic introners to matched
    non-focal background windows in the Group 1 4D all-sites VCF.
    """
    input:
        vcf = config["paths"]["fourfold_vcf"],
        targets = SWEEP_DIR / "group1_polymorphic_introner_targets.tsv"
    output:
        introner_windows = SWEEP_DIR / "introner_windows.tsv",
        background_windows = SWEEP_DIR / "background_windows.tsv",
        pvalues = SWEEP_DIR / "introner_sweep_empirical_pvalues.tsv",
        plot_png = FIGURES_DIR / "snp_popgen" / "group1_introner_sweep_summary.png",
        plot_pdf = FIGURES_DIR / "snp_popgen" / "group1_introner_sweep_summary.pdf"
    params:
        samples = ",".join(GROUP1_SAMPLES),
        window_sizes = ",".join(map(str, SWEEP_WINDOW_SIZES)),
        min_callable_sites = SWEEP_MIN_CALLABLE_SITES,
        backgrounds_per_focal = SWEEP_BACKGROUNDS_PER_FOCAL
    log:
        SELECTION_LOG_DIR / "calculate_group1_introner_sweep_windows.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {SWEEP_DIR} {FIGURES_DIR}/snp_popgen {SELECTION_LOG_DIR}

        python {PROJECT_ROOT}/scripts/popgen/selection_analysis/calculate_group1_introner_sweep_windows.py \
            --vcf {input.vcf} \
            --targets {input.targets} \
            --samples {params.samples} \
            --window-sizes {params.window_sizes} \
            --min-callable-sites {params.min_callable_sites} \
            --backgrounds-per-focal {params.backgrounds_per_focal} \
            --introner-windows {output.introner_windows} \
            --background-windows {output.background_windows} \
            --pvalues {output.pvalues} \
            --plot-png {output.plot_png} \
            --plot-pdf {output.plot_pdf} \
            > {log} 2>&1
        """


# ============================================================
# FIXED SNP SWEEP CONTROL
# ============================================================

rule build_fixed_snp_sweep_targets:
    """
    Build 4D SNP targets fixed between Group 1 and Group 2.

    These provide a same-timescale control for fixed introner loci: each target
    is a biallelic 4D SNP where all Group 1 samples carry one allele and all
    Group 2 samples carry the other allele.
    """
    input:
        vcf = config["paths"]["fourfold_all_samples_vcf"]
    output:
        tsv = FIXED_SNP_SWEEP_DIR / "fixed_snp_targets.tsv",
        bed = FIXED_SNP_SWEEP_DIR / "fixed_snp_targets.bed"
    params:
        group1 = ",".join(GROUP1_SAMPLES),
        group2 = ",".join(GROUP2_SAMPLES),
        min_spacing = SWEEP_FIXED_SNP_MIN_SPACING,
        mating_contig = lambda wildcards: f"{REFERENCE}#0#{config['mating_type_region']['scaffold']}",
        mating_start = config["mating_type_region"]["start"],
        mating_end = config["mating_type_region"]["end"]
    log:
        SELECTION_LOG_DIR / "build_fixed_snp_sweep_targets.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {FIXED_SNP_SWEEP_DIR} {SELECTION_LOG_DIR}

        python {PROJECT_ROOT}/scripts/popgen/selection_analysis/build_fixed_snp_sweep_targets.py \
            --vcf {input.vcf} \
            --output-tsv {output.tsv} \
            --output-bed {output.bed} \
            --group1-samples {params.group1} \
            --group2-samples {params.group2} \
            --min-spacing {params.min_spacing} \
            --mating-contig {params.mating_contig} \
            --mating-start {params.mating_start} \
            --mating-end {params.mating_end} \
            > {log} 2>&1
        """


rule calculate_fixed_snp_sweep_windows:
    """
    Compare sweep statistics around fixed 4D SNPs to matched non-focal
    background windows in the Group 1 4D all-sites VCF.
    """
    input:
        vcf = config["paths"]["fourfold_vcf"],
        targets = FIXED_SNP_SWEEP_DIR / "fixed_snp_targets.tsv"
    output:
        fixed_snp_windows = FIXED_SNP_SWEEP_DIR / "fixed_snp_windows.tsv",
        background_windows = FIXED_SNP_SWEEP_DIR / "fixed_snp_background_windows.tsv",
        pvalues = FIXED_SNP_SWEEP_DIR / "fixed_snp_sweep_empirical_pvalues.tsv",
        plot_png = FIGURES_DIR / "snp_popgen" / "fixed_snp_sweep_summary.png",
        plot_pdf = FIGURES_DIR / "snp_popgen" / "fixed_snp_sweep_summary.pdf"
    params:
        samples = ",".join(GROUP1_SAMPLES),
        window_sizes = ",".join(map(str, SWEEP_WINDOW_SIZES)),
        min_callable_sites = SWEEP_MIN_CALLABLE_SITES,
        backgrounds_per_focal = SWEEP_BACKGROUNDS_PER_FOCAL
    log:
        SELECTION_LOG_DIR / "calculate_fixed_snp_sweep_windows.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {FIXED_SNP_SWEEP_DIR} {FIGURES_DIR}/snp_popgen {SELECTION_LOG_DIR}

        python {PROJECT_ROOT}/scripts/popgen/selection_analysis/calculate_group1_introner_sweep_windows.py \
            --vcf {input.vcf} \
            --targets {input.targets} \
            --samples {params.samples} \
            --window-sizes {params.window_sizes} \
            --min-callable-sites {params.min_callable_sites} \
            --backgrounds-per-focal {params.backgrounds_per_focal} \
            --introner-windows {output.fixed_snp_windows} \
            --background-windows {output.background_windows} \
            --pvalues {output.pvalues} \
            --plot-png {output.plot_png} \
            --plot-pdf {output.plot_pdf} \
            --focal-label "fixed SNP" \
            > {log} 2>&1
        """


# ============================================================
# SWEEPFINDER2 GENOME-WIDE CLR SCAN
# ============================================================

rule build_sweepfinder_input:
    """
    Build SweepFinder2 allele-frequency and uniform grid files from the Group 1
    4D all-sites VCF (folded polymorphic sites, CCMP1545 reference frame).
    """
    input:
        vcf = config["paths"]["fourfold_vcf"],
        fai = CCMP1545_GRAPH_FAI
    output:
        combined_freq = SF2_DIR / "combined.freq",
        contigs = SF2_DIR / "contigs.tsv",
        freq_dir = directory(SF2_DIR / "freq"),
        grid_dir = directory(SF2_DIR / "grid")
    params:
        samples = ",".join(GROUP1_SAMPLES),
        grid_spacing = SF2_GRID_SPACING,
        min_called_haplotypes = SF2_MIN_CALLED_HAPLOTYPES,
        mating_contig = lambda wildcards: f"{REFERENCE}#0#{config['mating_type_region']['scaffold']}",
        mating_start = config["mating_type_region"]["start"],
        mating_end = config["mating_type_region"]["end"],
        exclude_contigs = lambda wildcards: f"{REFERENCE}#0#{config['mating_type_region']['scaffold']}"
    log:
        SELECTION_LOG_DIR / "build_sweepfinder_input.log"
    conda: "../envs/sweepfinder.yaml"
    shell:
        """
        mkdir -p {SF2_DIR}/freq {SF2_DIR}/grid {SELECTION_LOG_DIR}

        python {PROJECT_ROOT}/scripts/popgen/selection_analysis/build_sweepfinder_input.py \
            --vcf {input.vcf} \
            --fai {input.fai} \
            --samples {params.samples} \
            --output-dir {SF2_DIR} \
            --combined-freq {output.combined_freq} \
            --contigs-tsv {output.contigs} \
            --grid-spacing {params.grid_spacing} \
            --min-called-haplotypes {params.min_called_haplotypes} \
            --mating-contig {params.mating_contig} \
            --mating-start {params.mating_start} \
            --mating-end {params.mating_end} \
            --exclude-contigs {params.exclude_contigs} \
            > {log} 2>&1
        """


rule compute_sweepfinder_spectrum:
    """
    Compute the genome-wide empirical folded SFS for SweepFinder2 (-f).
    """
    input:
        combined_freq = SF2_DIR / "combined.freq"
    output:
        spectrum = SF2_DIR / "spectrum.out"
    log:
        SELECTION_LOG_DIR / "compute_sweepfinder_spectrum.log"
    conda: "../envs/sweepfinder.yaml"
    shell:
        """
        SweepFinder2 -f {input.combined_freq} {output.spectrum} > {log} 2>&1
        """


rule run_sweepfinder_scan:
    """
    Run SweepFinder2 CLR scan (-lu) on each contig using pre-computed spectrum.
    """
    input:
        contigs = SF2_DIR / "contigs.tsv",
        spectrum = SF2_DIR / "spectrum.out",
        freq_dir = SF2_DIR / "freq",
        grid_dir = SF2_DIR / "grid"
    output:
        marker = SF2_DIR / "scan_complete.txt",
        clr_dir = directory(SF2_DIR / "clr")
    log:
        SELECTION_LOG_DIR / "run_sweepfinder_scan.log"
    conda: "../envs/sweepfinder.yaml"
    shell:
        """
        mkdir -p {output.clr_dir}
        : > {log}

        tail -n +2 {input.contigs} | while IFS=$'\\t' read -r contig length n_sites n_grid freq_file grid_file; do
            safe=$(basename "$freq_file" .freq)
            echo "Scanning $contig ($n_grid grid points, $n_sites polymorphic sites)" >> {log}
            SweepFinder2 -lu "$grid_file" "$freq_file" {input.spectrum} {output.clr_dir}/${{safe}}.clr >> {log} 2>&1
        done

        echo "done" > {output.marker}
        """


rule call_sweep_regions:
    """
    Merge per-contig CLR outputs, call sweep-like regions, and plot Manhattan.
    """
    input:
        contigs = SF2_DIR / "contigs.tsv",
        marker = SF2_DIR / "scan_complete.txt",
        clr_dir = SF2_DIR / "clr"
    output:
        clr_tsv = SF2_DIR / "sweepfinder_clr.tsv",
        regions_bed = SF2_DIR / "sweep_regions.bed",
        summary = SF2_DIR / "sweep_calling_summary.tsv",
        plot_png = FIGURES_DIR / "snp_popgen" / "sweepfinder_manhattan.png",
        plot_pdf = FIGURES_DIR / "snp_popgen" / "sweepfinder_manhattan.pdf"
    params:
        clr_percentile = SF2_CLR_PERCENTILE,
        region_merge_gap = SF2_REGION_MERGE_GAP
    log:
        SELECTION_LOG_DIR / "call_sweep_regions.log"
    conda: "../envs/sweepfinder.yaml"
    shell:
        """
        mkdir -p {FIGURES_DIR}/snp_popgen {SELECTION_LOG_DIR}

        python {PROJECT_ROOT}/scripts/popgen/selection_analysis/call_sweep_regions.py \
            --clr-dir {input.clr_dir} \
            --contigs-tsv {input.contigs} \
            --output-clr {output.clr_tsv} \
            --output-regions {output.regions_bed} \
            --output-summary {output.summary} \
            --plot-png {output.plot_png} \
            --plot-pdf {output.plot_pdf} \
            --clr-percentile {params.clr_percentile} \
            --region-merge-gap {params.region_merge_gap} \
            > {log} 2>&1
        """


rule test_introner_sweep_enrichment:
    """
    Test fixed vs polymorphic Group 1 introner enrichment/depletion in
    SweepFinder2 sweep regions vs non-introner introns and a permutation null.
    """
    input:
        clr_tsv = SF2_DIR / "sweepfinder_clr.tsv",
        regions_bed = SF2_DIR / "sweep_regions.bed",
        contigs = SF2_DIR / "contigs.tsv",
        introner_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        non_introner_matrix = RESULTS / "evolution" / "non_introner_introns" / "intron_genotype_matrix.tsv"
    output:
        region_overlap = SF2_DIR / "introner_sweep_region_overlap.tsv",
        perlocus = SF2_DIR / "introner_sweep_perlocus.tsv",
        summary = SF2_DIR / "introner_sweep_enrichment_summary.tsv",
        enrichment_png = FIGURES_DIR / "snp_popgen" / "introner_sweep_enrichment.png",
        enrichment_pdf = FIGURES_DIR / "snp_popgen" / "introner_sweep_enrichment.pdf",
        clr_png = FIGURES_DIR / "snp_popgen" / "introner_sweep_clr_comparison.png",
        clr_pdf = FIGURES_DIR / "snp_popgen" / "introner_sweep_clr_comparison.pdf"
    params:
        reference_sample = REFERENCE,
        group1 = ",".join(GROUP1_SAMPLES),
        group2 = ",".join(GROUP2_SAMPLES),
        accepted_within_statuses = ",".join(SWEEP_ACCEPTED_WITHIN_STATUSES),
        flank_size = SF2_FLANK_SIZE,
        permutations = SF2_PERMUTATIONS,
        seed = SF2_SEED,
        mating_contig = lambda wildcards: f"{REFERENCE}#0#{config['mating_type_region']['scaffold']}",
        mating_start = config["mating_type_region"]["start"],
        mating_end = config["mating_type_region"]["end"]
    log:
        SELECTION_LOG_DIR / "test_introner_sweep_enrichment.log"
    conda: "../envs/sweepfinder.yaml"
    shell:
        """
        mkdir -p {FIGURES_DIR}/snp_popgen {SELECTION_LOG_DIR}

        python {PROJECT_ROOT}/scripts/popgen/selection_analysis/test_introner_sweep_enrichment.py \
            --clr-tsv {input.clr_tsv} \
            --sweep-regions {input.regions_bed} \
            --contigs-tsv {input.contigs} \
            --introner-matrix {input.introner_matrix} \
            --non-introner-matrix {input.non_introner_matrix} \
            --group1-samples {params.group1} \
            --group2-samples {params.group2} \
            --accepted-within-statuses {params.accepted_within_statuses} \
            --reference-sample {params.reference_sample} \
            --flank-size {params.flank_size} \
            --mating-contig {params.mating_contig} \
            --mating-start {params.mating_start} \
            --mating-end {params.mating_end} \
            --permutations {params.permutations} \
            --seed {params.seed} \
            --region-overlap-tsv {output.region_overlap} \
            --perlocus-tsv {output.perlocus} \
            --summary-tsv {output.summary} \
            --plot-enrichment-png {output.enrichment_png} \
            --plot-enrichment-pdf {output.enrichment_pdf} \
            --plot-clr-png {output.clr_png} \
            --plot-clr-pdf {output.clr_pdf} \
            > {log} 2>&1
        """


# ============================================================
# PRINCIPAL COMPONENT ANALYSIS
# ============================================================

rule vcf_to_plink:
    """
    Convert VCF to PLINK format for PCA.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz"
    output:
        bed = PCA_DIR / "mpusilla.snps.bed",
        bim = PCA_DIR / "mpusilla.snps.bim",
        fam = PCA_DIR / "mpusilla.snps.fam"
    params:
        prefix = PCA_DIR / "mpusilla.snps"
    log:
        SELECTION_LOG_DIR / "vcf_to_plink.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {PCA_DIR}

        plink --vcf {input.vcf} \
            --make-bed \
            --out {params.prefix} \
            --allow-extra-chr \
            --double-id \
            2> {log}
        """


rule ld_pruning:
    """
    LD pruning for PCA to remove correlated SNPs.

    Uses sliding window approach:
    - Window size: 50 SNPs
    - Step size: 5 SNPs
    - r² threshold: 0.2
    """
    input:
        bed = PCA_DIR / "mpusilla.snps.bed",
        bim = PCA_DIR / "mpusilla.snps.bim",
        fam = PCA_DIR / "mpusilla.snps.fam"
    output:
        prune_in = PCA_DIR / "mpusilla.snps.prune.in",
        prune_out = PCA_DIR / "mpusilla.snps.prune.out"
    params:
        prefix = PCA_DIR / "mpusilla.snps",
        window = 50,
        step = 5,
        r2 = 0.2
    log:
        SELECTION_LOG_DIR / "ld_pruning.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        plink --bfile {params.prefix} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --out {params.prefix} \
            --allow-extra-chr \
            2> {log}
        """


rule run_pca:
    """
    Run PCA on LD-pruned SNPs.
    """
    input:
        bed = PCA_DIR / "mpusilla.snps.bed",
        prune_in = PCA_DIR / "mpusilla.snps.prune.in"
    output:
        eigenval = PCA_DIR / "mpusilla.snps.eigenval",
        eigenvec = PCA_DIR / "mpusilla.snps.eigenvec"
    params:
        prefix = PCA_DIR / "mpusilla.snps"
    log:
        SELECTION_LOG_DIR / "pca.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        plink --bfile {params.prefix} \
            --extract {input.prune_in} \
            --pca 10 \
            --out {params.prefix} \
            --allow-extra-chr \
            2> {log}
        """


rule plot_pca:
    """
    Generate PCA visualization plot.

    Colors samples by group:
    - Group1 (intronerful): Blue
    - Group2 (intronerless): Red
    """
    input:
        eigenval = PCA_DIR / "mpusilla.snps.eigenval",
        eigenvec = PCA_DIR / "mpusilla.snps.eigenvec"
    output:
        pdf = FIGURES_DIR / "snp_popgen" / "pca.pdf",
        png = FIGURES_DIR / "snp_popgen" / "pca.png"
    params:
        group1_str = ",".join(GROUP1_SAMPLES),
        group2_str = ",".join(GROUP2_SAMPLES)
    log:
        SELECTION_LOG_DIR / "plot_pca.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/pca/plot_PCA.py \
            --eigenval {input.eigenval} \
            --eigenvec {input.eigenvec} \
            --output_pdf {output.pdf} \
            --output_png {output.png} \
            --group1 {params.group1_str} \
            --group2 {params.group2_str} \
            2> {log}
        """


# ============================================================
# LINKAGE DISEQUILIBRIUM ANALYSIS
# ============================================================

rule compute_ld:
    """
    Compute pairwise linkage disequilibrium (r²) between SNPs.

    Calculates LD within 100kb windows.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.group1.vcf.gz"
    output:
        ld = LD_DIR / "ld_pairwise.tsv.gz"
    params:
        max_distance = 100000  # 100kb
    log:
        SELECTION_LOG_DIR / "compute_ld.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {LD_DIR}

        python {PROJECT_ROOT}/scripts/popgen/ld/compute_ld.py \
            --vcf {input.vcf} \
            --output {output.ld} \
            --max_distance {params.max_distance} \
            2> {log}
        """


rule summarize_ld:
    """
    Summarize LD by distance bins.

    Groups SNP pairs into distance bins and calculates
    mean r² per bin for LD decay visualization.
    """
    input:
        ld = LD_DIR / "ld_pairwise.tsv.gz"
    output:
        summary = LD_DIR / "ld_summary_by_distance.tsv"
    params:
        bin_size = 10000  # 10kb bins
    log:
        SELECTION_LOG_DIR / "summarize_ld.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/ld/summarize_ld.py \
            {input.ld} \
            {output.summary} \
            --bin_size {params.bin_size} \
            2> {log}
        """


rule plot_ld_decay:
    """
    Plot LD decay with genomic distance.
    """
    input:
        summary = LD_DIR / "ld_summary_by_distance.tsv"
    output:
        pdf = FIGURES_DIR / "snp_popgen" / "ld_decay.pdf",
        png = FIGURES_DIR / "snp_popgen" / "ld_decay.png"
    log:
        SELECTION_LOG_DIR / "plot_ld.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/ld/plot_ld.py \
            {input.summary} \
            {output.pdf} \
            2> {log}

        python {PROJECT_ROOT}/scripts/popgen/ld/plot_ld.py \
            {input.summary} \
            {output.png} \
            2>> {log}
        """


# ============================================================
# RECOMBINATION ANALYSIS
# ============================================================

rule compare_recombination_introners_all:
    """
    Gene-centric recombination comparison: all introners vs non-introner windows.
    Produces gene_exonic_introners_all_5kb_updated_{summary,windows}.tsv + pdf.
    """
    input:
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        introner_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        summary_tsv = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated_summary.tsv",
        windows_tsv = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated_windows.tsv",
        pdf = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated.pdf",
        png = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated.png"
    params:
        pyrho_dir = PYRHO_CCMP1545,
        sample_name = "CCMP1545",
        introner_type = "all",
        window_size = 5000,
        merge_distance = 10000,
        output_prefix = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated"
    log:
        SELECTION_LOG_DIR / "recombination_introners_all.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {RECOMB_DIR}

        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/compare_recombination_gene_exonic_introners.py \
            --pyrho_dir {params.pyrho_dir} \
            --gtf_file {input.gtf} \
            --introner_matrix {input.introner_matrix} \
            --sample_name {params.sample_name} \
            --introner_type {params.introner_type} \
            --window_size {params.window_size} \
            --merge_distance {params.merge_distance} \
            --output_prefix {params.output_prefix} \
            2> {log}
        """


rule compare_recombination_introners_polymorphic:
    """
    Gene-centric recombination comparison: polymorphic introners only.
    Produces gene_exonic_polymorphic_introners_5kb_updated_{summary,windows}.tsv + pdf.
    """
    input:
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        introner_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        summary_tsv = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated_summary.tsv",
        windows_tsv = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated_windows.tsv",
        pdf = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated.pdf",
        png = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated.png"
    params:
        pyrho_dir = PYRHO_CCMP1545,
        sample_name = "CCMP1545",
        introner_type = "polymorphic",
        window_size = 5000,
        merge_distance = 10000,
        output_prefix = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated"
    log:
        SELECTION_LOG_DIR / "recombination_introners_polymorphic.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {RECOMB_DIR}

        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/compare_recombination_gene_exonic_introners.py \
            --pyrho_dir {params.pyrho_dir} \
            --gtf_file {input.gtf} \
            --introner_matrix {input.introner_matrix} \
            --sample_name {params.sample_name} \
            --introner_type {params.introner_type} \
            --window_size {params.window_size} \
            --merge_distance {params.merge_distance} \
            --output_prefix {params.output_prefix} \
            2> {log}
        """


rule compare_recombination_frequency_based:
    """
    Frequency-based recombination comparison: recent gain vs recent loss vs non-introner.
    Produces gene_exonic_frequency_based_5kb_updated_{summary,windows}.tsv + pdf.
    """
    input:
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        introner_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        summary_tsv = RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated_summary.tsv",
        windows_tsv = RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated_windows.tsv",
        pdf = RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated.pdf",
        png = RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated.png"
    params:
        pyrho_dir = PYRHO_CCMP1545,
        exclude_samples = ",".join(GROUP2_SAMPLES),
        window_size = 5000,
        merge_distance = 10000,
        output_prefix = RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated"
    log:
        SELECTION_LOG_DIR / "recombination_frequency_based.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {RECOMB_DIR}

        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/compare_recombination_gene_exonic_frequency_based.py \
            --pyrho_dir {params.pyrho_dir} \
            --gtf_file {input.gtf} \
            --introner_matrix {input.introner_matrix} \
            --exclude_samples {params.exclude_samples} \
            --window_size {params.window_size} \
            --merge_distance {params.merge_distance} \
            --output_prefix {params.output_prefix} \
            2> {log}
        """


rule compare_recombination_group2:
    """
    Recombination at positions where Group 2 (RCC1749/RCC3052) has introners,
    measured in RCC1749's recombination landscape.

    Uses the RCC1749 GTF (contigs named RCC1749#0#intronerless_contig_X) and
    RCC1749 pyrho map together — both are in RCC1749 frame and align with the
    G2-sample contig names in the genotype matrix.
    """
    input:
        gtf = ANNOTATIONS_DIR / "RCC1749.gtf",
        introner_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        summary_tsv = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated_summary.tsv",
        windows_tsv = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated_windows.tsv",
        pdf = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated.pdf",
        png = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated.png"
    params:
        pyrho_dir = PYRHO_RCC1749,
        window_size = 5000,
        merge_distance = 10000,
        output_prefix = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated"
    log:
        SELECTION_LOG_DIR / "recombination_group2.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {RECOMB_DIR}

        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/compare_recombination_gene_exonic_group2_introners.py \
            --pyrho_dir {params.pyrho_dir} \
            --gtf_file {input.gtf} \
            --introner_matrix {input.introner_matrix} \
            --window_size {params.window_size} \
            --merge_distance {params.merge_distance} \
            --output_prefix {params.output_prefix} \
            2> {log}
        """


rule plot_recombination_boxplots:
    """
    Generate boxplots comparing recombination rates by introner category:
    non-introner, all introners, polymorphic, recent gain, recent loss, plus Group 2.
    """
    input:
        all_windows = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated_windows.tsv",
        poly_windows = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated_windows.tsv",
        freq_windows = RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated_windows.tsv",
        group2_windows = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated_windows.tsv"
    output:
        pdf = FIGURES_DIR / "snp_popgen" / "recombination_boxplots.pdf",
        png = FIGURES_DIR / "snp_popgen" / "recombination_boxplots.png"
    params:
        base_dir = RECOMB_DIR,
        output_prefix = FIGURES_DIR / "snp_popgen" / "recombination_boxplots"
    log:
        SELECTION_LOG_DIR / "plot_recombination.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {FIGURES_DIR}/snp_popgen

        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/plot_recombination_boxplots.py \
            --base_dir {params.base_dir} \
            --output {params.output_prefix} \
            2> {log}
        """


rule plot_recombination_violin_categories:
    """
    Generate violin plots for selected pyrho recombination categories:
    Population 1 exonic background, all introners, polymorphic introners,
    Population 2 exonic background, and Population 2 all introners.
    """
    input:
        group1_all_windows = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated_windows.tsv",
        group1_polymorphic_windows = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated_windows.tsv",
        group2_windows = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated_windows.tsv",
        color_guide = PROJECT_ROOT / "master_figure_color_guide.tsv",
        script = PROJECT_ROOT / "scripts" / "popgen" / "recombination_analysis" / "plot_recombination_violin_categories.py"
    output:
        pdf = FIGURES_DIR / "snp_popgen" / "recombination_violin_categories.pdf",
        png = FIGURES_DIR / "snp_popgen" / "recombination_violin_categories.png",
        tsv = RECOMB_DIR / "recombination_violin_categories.tsv"
    log:
        SELECTION_LOG_DIR / "plot_recombination_violin_categories.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {FIGURES_DIR}/snp_popgen {RECOMB_DIR}

        python {input.script} \
            --group1-all-windows {input.group1_all_windows} \
            --group1-polymorphic-windows {input.group1_polymorphic_windows} \
            --group2-windows {input.group2_windows} \
            --output-pdf {output.pdf} \
            --output-png {output.png} \
            --output-tsv {output.tsv} \
            --color-guide {input.color_guide} \
            2> {log}
        """


# ============================================================
# TARGET RULES
# ============================================================

rule selection_complete:
    """
    Target: Complete selection and population structure analysis.
    """
    input:
        # Basic popgen
        SELECTION_DIR / "tajimas_d.tsv",
        # PCA
        PCA_DIR / "mpusilla.snps.eigenvec",
        FIGURES_DIR / "snp_popgen" / "pca.pdf",
        # LD
        LD_DIR / "ld_summary_by_distance.tsv",
        FIGURES_DIR / "snp_popgen" / "ld_decay.pdf",
        # Recombination
        RECOMB_DIR / "gene_exonic_introners_all_5kb_updated_summary.tsv",
        RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated_summary.tsv",
        RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated_summary.tsv",
        RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated_summary.tsv",
        FIGURES_DIR / "snp_popgen" / "recombination_boxplots.pdf",
        FIGURES_DIR / "snp_popgen" / "recombination_violin_categories.pdf",
        # SFS
        SFS_DIR / "afs_by_class.tsv",
        FIGURES_DIR / "snp_popgen" / "afs_by_class.pdf",
        # SweepFinder2
        SF2_DIR / "sweepfinder_clr.tsv",
        SF2_DIR / "sweep_regions.bed",
        SF2_DIR / "introner_sweep_enrichment_summary.tsv",
        FIGURES_DIR / "snp_popgen" / "sweepfinder_manhattan.png",
        FIGURES_DIR / "snp_popgen" / "introner_sweep_enrichment.png",


rule basic_popgen_only:
    """
    Target: Basic population genetics statistics only.
    """
    input:
        SELECTION_DIR / "tajimas_d.tsv",
        SFS_DIR / "afs_by_class.tsv",
        FIGURES_DIR / "snp_popgen" / "afs_by_class.pdf"


rule pca_only:
    """
    Target: PCA analysis only.
    """
    input:
        PCA_DIR / "mpusilla.snps.eigenvec",
        FIGURES_DIR / "snp_popgen" / "pca.pdf"


rule ld_only:
    """
    Target: LD analysis only.
    """
    input:
        LD_DIR / "ld_summary_by_distance.tsv",
        FIGURES_DIR / "snp_popgen" / "ld_decay.pdf"


rule recombination_only:
    """
    Target: Recombination analysis only.
    """
    input:
        RECOMB_DIR / "gene_exonic_introners_all_5kb_updated_summary.tsv",
        RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated_summary.tsv",
        RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated_summary.tsv",
        RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated_summary.tsv",
        FIGURES_DIR / "snp_popgen" / "recombination_boxplots.pdf",
        FIGURES_DIR / "snp_popgen" / "recombination_violin_categories.pdf"


rule group1_introner_sweeps:
    """
    Target: Group 1 polymorphic introner selective sweep screen only.
    """
    input:
        SWEEP_DIR / "group1_polymorphic_introner_targets.tsv",
        SWEEP_DIR / "group1_polymorphic_introner_targets.bed",
        SWEEP_DIR / "introner_windows.tsv",
        SWEEP_DIR / "background_windows.tsv",
        SWEEP_DIR / "introner_sweep_empirical_pvalues.tsv",
        FIGURES_DIR / "snp_popgen" / "group1_introner_sweep_summary.png",
        FIGURES_DIR / "snp_popgen" / "group1_introner_sweep_summary.pdf"


rule fixed_snp_sweep_control:
    """
    Target: Fixed 4D SNP sweep-control screen only.
    """
    input:
        FIXED_SNP_SWEEP_DIR / "fixed_snp_targets.tsv",
        FIXED_SNP_SWEEP_DIR / "fixed_snp_targets.bed",
        FIXED_SNP_SWEEP_DIR / "fixed_snp_windows.tsv",
        FIXED_SNP_SWEEP_DIR / "fixed_snp_background_windows.tsv",
        FIXED_SNP_SWEEP_DIR / "fixed_snp_sweep_empirical_pvalues.tsv",
        FIGURES_DIR / "snp_popgen" / "fixed_snp_sweep_summary.png",
        FIGURES_DIR / "snp_popgen" / "fixed_snp_sweep_summary.pdf"


rule sweepfinder_sweeps:
    """
    Target: SweepFinder2 genome-wide CLR scan and introner enrichment test.
    """
    input:
        SF2_DIR / "combined.freq",
        SF2_DIR / "contigs.tsv",
        SF2_DIR / "spectrum.out",
        SF2_DIR / "scan_complete.txt",
        SF2_DIR / "sweepfinder_clr.tsv",
        SF2_DIR / "sweep_regions.bed",
        SF2_DIR / "sweep_calling_summary.tsv",
        SF2_DIR / "introner_sweep_region_overlap.tsv",
        SF2_DIR / "introner_sweep_perlocus.tsv",
        SF2_DIR / "introner_sweep_enrichment_summary.tsv",
        FIGURES_DIR / "snp_popgen" / "sweepfinder_manhattan.png",
        FIGURES_DIR / "snp_popgen" / "sweepfinder_manhattan.pdf",
        FIGURES_DIR / "snp_popgen" / "introner_sweep_enrichment.png",
        FIGURES_DIR / "snp_popgen" / "introner_sweep_enrichment.pdf",
        FIGURES_DIR / "snp_popgen" / "introner_sweep_clr_comparison.png",
        FIGURES_DIR / "snp_popgen" / "introner_sweep_clr_comparison.pdf"
