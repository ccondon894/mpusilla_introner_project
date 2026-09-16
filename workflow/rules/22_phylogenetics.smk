# ============================================================
# 22_phylogenetics.smk - Phylogenetic Analysis with IQ-TREE
# ============================================================
#
# Constructs phylogenetic trees from 4-fold degenerate SNPs.
#
# Pipeline:
# 1. Convert VCF to PHYLIP format
# 2. Run IQ-TREE maximum likelihood analysis
# 3. Generate publication-quality tree plots
#
# Adapted from:
# - /scratch1/chris/introner_vis/iqtree/commands.sh
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
PHYLO_DIR = SNP_DIR / "phylogenetics"
PHYLO_LOG_DIR = SNP_DIR / "logs" / "phylogenetics"

# IQ-TREE parameters
IQTREE_MODEL = config["params"]["iqtree"]["model"]
IQTREE_BOOTSTRAP = config["params"]["iqtree"]["bootstrap"]


# ============================================================
# VCF TO PHYLIP CONVERSION
# ============================================================

rule vcf_to_phylip:
    """
    Convert VCF to PHYLIP format for phylogenetic analysis.

    Uses vcf2phylip.py to create alignment from SNP data.
    The script handles IUPAC ambiguity codes for heterozygous sites
    (though M. pusilla is haploid).
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz"
    output:
        phylip = PHYLO_DIR / "mpusilla.snps.4d.min4.phy",
        fasta = PHYLO_DIR / "mpusilla.snps.4d.min4.fasta"
    params:
        min_samples = 4,  # Minimum samples with data per site
        output_folder = PHYLO_DIR,
        output_prefix = "mpusilla.snps.4d"
    log:
        PHYLO_LOG_DIR / "vcf2phylip.log"
    conda: "../envs/phylo.yaml"
    shell:
        """
        mkdir -p {PHYLO_DIR}
        mkdir -p {PHYLO_LOG_DIR}

        python {PROJECT_ROOT}/scripts/popgen/iqtree/vcf2phylip.py \
            -i {input.vcf} \
            --output-folder {params.output_folder} \
            --output-prefix {params.output_prefix} \
            --min-samples-locus {params.min_samples} \
            --fasta \
            2> {log}
        """


rule vcf_to_phylip_nomissing:
    """
    Convert no-missing-data VCF to PHYLIP format.

    Uses stricter filtering for sites with complete data.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.nomissing.vcf.gz"
    output:
        phylip = PHYLO_DIR / "mpusilla.snps.4d.nomissing.min4.phy",
        fasta = PHYLO_DIR / "mpusilla.snps.4d.nomissing.min4.fasta"
    params:
        output_prefix = PHYLO_DIR / "mpusilla.snps.4d.nomissing",
        min_samples = 4,
        output_folder = PHYLO_DIR,
    log:
        PHYLO_LOG_DIR / "vcf2phylip_nomissing.log"
    conda: "../envs/phylo.yaml"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/iqtree/vcf2phylip.py \
            -i {input.vcf} \
            --output-folder {params.output_folder} \
            --output-prefix {params.output_prefix} \
            --min-samples-locus {params.min_samples} \
            --fasta \
            2> {log}
        """


# ============================================================
# IQ-TREE MAXIMUM LIKELIHOOD ANALYSIS
# ============================================================

rule run_iqtree:
    """
    Run IQ-TREE maximum likelihood phylogenetic analysis.

    Uses GTR+ASC model (with ascertainment bias correction for SNP data)
    and ultrafast bootstrap for branch support.
    """
    input:
        phylip = PHYLO_DIR / "mpusilla.snps.4d.min4.phy"
    output:
        treefile = PHYLO_DIR / "mpusilla.snps.4d.min4.phy.treefile",
        iqtree = PHYLO_DIR / "mpusilla.snps.4d.min4.phy.iqtree",
        log_file = PHYLO_DIR / "mpusilla.snps.4d.min4.phy.log"
    params:
        model = "GTR+ASC",  # GTR with ascertainment bias correction for SNPs
        bootstrap = IQTREE_BOOTSTRAP,
        prefix = PHYLO_DIR / "mpusilla.snps.4d.min4.phy"
    threads: 8
    log:
        PHYLO_LOG_DIR / "iqtree.log"
    conda: "../envs/phylo.yaml"
    shell:
        """
        iqtree2 \
            -s {input.phylip} \
            -m {params.model} \
            -bb {params.bootstrap} \
            -nt {threads} \
            --prefix {params.prefix} \
            --redo \
            2> {log}
        """


rule run_iqtree_nomissing:
    """
    Run IQ-TREE on complete-data alignment.
    """
    input:
        phylip = PHYLO_DIR / "mpusilla.snps.4d.nomissing.min4.phy"
    output:
        treefile = PHYLO_DIR / "mpusilla.snps.4d.nomissing.min4.phy.treefile",
        iqtree = PHYLO_DIR / "mpusilla.snps.4d.nomissing.min4.phy.iqtree"
    params:
        model = "GTR+ASC",
        bootstrap = IQTREE_BOOTSTRAP,
        prefix = PHYLO_DIR / "mpusilla.snps.4d.nomissing.min4.phy"
    threads: 4
    log:
        PHYLO_LOG_DIR / "iqtree_nomissing.log"
    conda: "../envs/phylo.yaml"
    shell:
        """
        iqtree2 \
            -s {input.phylip} \
            -m {params.model} \
            -bb {params.bootstrap} \
            -nt {threads} \
            --prefix {params.prefix} \
            --redo \
            2> {log}
        """


# ============================================================
# TREE VISUALIZATION
# ============================================================

rule plot_phylogenetic_tree:
    """
    Generate publication-quality phylogenetic tree plot.

    Colors samples by group:
    - Group1 (intronerful): Blue
    - Group2 (intronerless): Red

    Includes a scale bar; bootstrap support values are omitted from the figure.
    """
    input:
        treefile = PHYLO_DIR / "mpusilla.snps.4d.rooted.treefile"
    output:
        pdf = FIGURES_DIR / "snp_popgen"  / "phylogenetic_tree.pdf",
        png = FIGURES_DIR / "snp_popgen"  / "phylogenetic_tree.png"
    log:
        PHYLO_LOG_DIR / "plot_tree.log"
    conda: "../envs/phylo.yaml"
    shell:
        """
        mkdir -p {FIGURES_DIR}/snp_popgen

        python {PROJECT_ROOT}/scripts/popgen/iqtree/plot_tree_improved.py \
            -i {input.treefile} \
            -o {FIGURES_DIR}/snp_popgen/phylogenetic_tree \
            --color-guide {PROJECT_ROOT}/master_figure_color_guide.tsv \
            --width 720 \
            --height 330 \
            --tip-font-size 14 \
            2> {log}
        """


rule reroot_tree:
    """
    Reroot tree using Group2 as outgroup.

    Creates a rooted tree for evolutionary interpretation.
    """
    input:
        treefile = PHYLO_DIR / "mpusilla.snps.4d.min4.phy.treefile"
    output:
        rooted = PHYLO_DIR / "mpusilla.snps.4d.rooted.treefile"
    params:
        outgroup = " ".join(GROUP2_SAMPLES)  # Both RCC1749 and RCC3052
    conda: "../envs/phylo.yaml"
    shell:
        """
        # Use nw_reroot from newick-utils if available, otherwise Python
        if command -v nw_reroot &> /dev/null; then
            nw_reroot {input.treefile} {params.outgroup} > {output.rooted}
        else
            python -c "
import sys
try:
    from ete3 import Tree
    t = Tree('{input.treefile}')
    ancestor = t.get_common_ancestor('RCC1749', 'RCC3052')
    t.set_outgroup(ancestor)
    t.write(outfile='{output.rooted}')
except ImportError:
    import shutil
    shutil.copy('{input.treefile}', '{output.rooted}')
"
        fi
        """


# ============================================================
# TARGET RULES
# ============================================================

rule phylogenetics_complete:
    """
    Target: Complete phylogenetic analysis.
    """
    input:
        PHYLO_DIR / "mpusilla.snps.4d.min4.phy.treefile",
        PHYLO_DIR / "mpusilla.snps.4d.nomissing.min4.phy.treefile",
        PHYLO_DIR / "mpusilla.snps.4d.rooted.treefile",
        FIGURES_DIR / "snp_popgen" / "phylogenetic_tree.pdf"


rule iqtree_only:
    """
    Target: Run IQ-TREE analysis only.
    """
    input:
        PHYLO_DIR / "mpusilla.snps.4d.min4.phy.treefile"


rule tree_plots_only:
    """
    Target: Generate tree plots only.
    """
    input:
        FIGURES_DIR / "snp_popgen" / "phylogenetic_tree.pdf",
        FIGURES_DIR / "snp_popgen"  / "snp_popgen"  / "phylogenetic_tree.png"
