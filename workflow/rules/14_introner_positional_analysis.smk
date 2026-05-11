# ============================================================================
# 14_introner_positional_analysis.smk
# Spearman correlation of introner density across genomic windows
# Inter-genome (CCMP1545 vs RCC1749) and intra-genome (family vs family)
#
# Source: /scratch1/chris/introner-group-comparison-analysis/
# ============================================================================

import os

# Strain identifiers for the inter-genome comparison
G1 = config["samples"]["reference"]["group1"]   # CCMP1545
G2 = config["samples"]["reference"]["group2"]   # RCC1749

# Output subdirectories
INTRONER_POS_DIR = RESULTS / "introner_positional"
INTRONER_POS_FIGURES = FIGURES_DIR / "introner_positional"

# Script directory
IP_SCRIPTS = str(Path(workflow.basedir).parent / "scripts" / "introner_positional")

# Config params
IP_PARAMS = config.get("params", {}).get("introner_positional", {})
WINDOW_SIZES_KB = IP_PARAMS.get("window_sizes_kb", [1, 5, 10, 20, 50])
IP_FAMILIES = IP_PARAMS.get("families", ['all', '1', '2', '3', '4', '14', '15'])
INTRA_FAMILY_PAIRS = IP_PARAMS.get("intra_genome_family_pairs", [[2, 4]])
INTRA_GENOMES = IP_PARAMS.get("intra_genome_genomes", ['CCMP1545', 'RCC1749'])
MIN_ALN_COV = IP_PARAMS.get("min_alignment_coverage", 0.5)
EXCLUDE_FIXED = IP_PARAMS.get("exclude_fixed", True)


def ip_family_label(family):
    """Convert family name to filename-safe label."""
    return family if family == 'all' else f'family_{family}'


# ============================================================================
# TARGET RULE
# ============================================================================

rule introner_positional_all:
    """
    Target: Introner positional correlation analysis

    Produces inter-genome and intra-genome Spearman correlation
    analyses of introner density across window scales.
    """
    input:
        INTRONER_POS_FIGURES / "inter_genome_window_scale_correlation.png"


# ============================================================================
# INTER-GENOME CORRELATION RULES
# ============================================================================

rule inter_genome_spearman_correlation:
    """Compute inter-genome Spearman correlation of introner density."""
    input:
        delta = GENOME_ALIGNMENT_DIR / "mummer" / f"{G1}_vs_{G2}.1delta",
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        g1_fai = ASSEMBLIES_DIR / f"{G1}.vg_paths.fa.fai",
        g2_fai = ASSEMBLIES_DIR / f"{G2}.vg_paths.fa.fai"
    output:
        json_out = INTRONER_POS_DIR / "inter_genome" / "spearman_{ip_family}_{window_kb}kb.json",
        plot_out = INTRONER_POS_FIGURES / "inter_genome" / "hexbin_{ip_family}_{window_kb}kb.png"
    params:
        family = lambda wc: wc.ip_family.replace('family_', '') if wc.ip_family != 'all' else 'all',
        min_cov = MIN_ALN_COV,
        exclude_fixed = "--exclude-fixed" if EXCLUDE_FIXED else ""
    log:
        INTRONER_POS_DIR / "logs" / "inter_genome_{ip_family}_{window_kb}kb.log"
    shell:
        """
        mkdir -p $(dirname {output.json_out}) $(dirname {output.plot_out}) $(dirname {log}) && \
        python3 {IP_SCRIPTS}/genome_wide_introner_correlation.py \
            --delta {input.delta} \
            --genotype-matrix {input.genotype_matrix} \
            --g1-fai {input.g1_fai} \
            --g2-fai {input.g2_fai} \
            --g1-name {G1} \
            --g2-name {G2} \
            --window-size $(({wildcards.window_kb} * 1000)) \
            --introner-family {params.family} \
            --min-alignment-coverage {params.min_cov} \
            {params.exclude_fixed} \
            --output {output.json_out} \
            --plot {output.plot_out} \
            2>&1 | tee {log}
        """


rule plot_inter_genome_window_scale:
    """Summary line plot of inter-genome correlation across window scales."""
    input:
        jsons = expand(
            INTRONER_POS_DIR / "inter_genome" / "spearman_{fl}_{wk}kb.json",
            fl=[ip_family_label(f) for f in IP_FAMILIES],
            wk=WINDOW_SIZES_KB
        )
    output:
        INTRONER_POS_FIGURES / "inter_genome_window_scale_correlation.png"
    params:
        result_dir = str(INTRONER_POS_DIR / "inter_genome"),
        families = ' '.join(str(f) for f in IP_FAMILIES),
        window_sizes = ' '.join(str(w) for w in WINDOW_SIZES_KB)
    log:
        INTRONER_POS_DIR / "logs" / "plot_inter_genome_window_scale.log"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log}) && \
        python3 {IP_SCRIPTS}/plot_window_scale_correlation.py \
            --result-dir {params.result_dir} \
            --families {params.families} \
            --window-sizes {params.window_sizes} \
            --output {output} \
            2>&1 | tee {log}
        """


# ============================================================================
# INTRA-GENOME CORRELATION RULES
# ============================================================================

rule intra_genome_spearman_correlation:
    """Compute intra-genome Spearman correlation between introner families."""
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        genome_fai = lambda wc: ASSEMBLIES_DIR / f"{wc.ip_genome}.vg_paths.fa.fai"
    output:
        json_out = INTRONER_POS_DIR / "intra_genome" / "spearman_{ip_genome}_{fam1}_vs_{fam2}_{window_kb}kb.json",
        plot_out = INTRONER_POS_FIGURES / "intra_genome" / "hexbin_{ip_genome}_{fam1}_vs_{fam2}_{window_kb}kb.png"
    params:
        exclude_fixed = "--exclude-fixed" if EXCLUDE_FIXED else ""
    log:
        INTRONER_POS_DIR / "logs" / "intra_genome_{ip_genome}_{fam1}_vs_{fam2}_{window_kb}kb.log"
    shell:
        """
        mkdir -p $(dirname {output.json_out}) $(dirname {output.plot_out}) $(dirname {log}) && \
        python3 {IP_SCRIPTS}/intra_genome_family_correlation.py \
            --genotype-matrix {input.genotype_matrix} \
            --genome-fai {input.genome_fai} \
            --genome-name {wildcards.ip_genome} \
            --family1 {wildcards.fam1} \
            --family2 {wildcards.fam2} \
            --window-size $(({wildcards.window_kb} * 1000)) \
            {params.exclude_fixed} \
            --output {output.json_out} \
            --plot {output.plot_out} \
            2>&1 | tee {log}
        """
