# ============================================================
# 40_go_enrichment.smk - GO Term Enrichment Analysis
# ============================================================
#
# Performs Gene Ontology enrichment analysis for introner genes:
# 1. Generate GO mappings from multiple databases (Pfam, KO, TAIR, PANTHER)
# 2. Enhance GO coverage by merging all sources
# 3. Phase 1: Gene classification by introner status
# 4. Phase 2: GO enrichment analysis (Fixed vs Polymorphic, Group1 vs Group2)
# 5. Phase 2b: Family-stratified enrichment
# 6. Visualization and reporting
#
# Adapted from:
# - /scratch1/chris/mpusilla_go_analysis/scripts/
#
# ============================================================

import os
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Output directories
GO_DIR = RESULTS / "go_enrichment"
GO_MAPPING_DIR = GO_DIR / "mappings"
GO_RESULTS_DIR = GO_DIR / "results"
GO_PLOTS_DIR = GO_DIR / "plots"
GO_LOG_DIR = GO_DIR / "logs"

# External database paths
INTERPRO2GO = config["paths"]["external_db"]["interpro2go"]
TAIR_ASSOC = config["paths"]["external_db"]["gene_association_tair"]
PANTHER_HMM = config["paths"]["external_db"]["panther_hmm"]

# Annotation file
ANNOTATION_INFO = config["paths"]["references"]["annotation_info"]


# ============================================================
# GO TERM MAPPING FROM MULTIPLE DATABASES
# ============================================================

rule pfam_to_go:
    """
    Map Pfam domain IDs to GO terms via InterPro.

    Uses InterPro2GO mapping file and queries InterPro API
    to link Pfam domains to GO terms.
    """
    input:
        interpro2go = INTERPRO2GO,
        annotation = ANNOTATION_INFO
    output:
        json = GO_MAPPING_DIR / "pfam_to_go.json"
    log:
        GO_LOG_DIR / "pfam_to_go.log"
    shell:
        """
        mkdir -p {GO_MAPPING_DIR}
        mkdir -p {GO_LOG_DIR}

        python {PROJECT_ROOT}/scripts/go_analysis/pfam_to_go.py \
            {input.interpro2go} \
            {input.annotation} \
            {output.json} \
            2> {log}
        """


rule ko_to_go:
    """
    Map KEGG Orthology (KO) IDs to GO terms via KEGG REST API.

    NOTE: Requires internet connection.
    """
    input:
        annotation = ANNOTATION_INFO
    output:
        json = GO_MAPPING_DIR / "ko_to_go.json"
    log:
        GO_LOG_DIR / "ko_to_go.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/go_analysis/ko_to_go.py \
            {input.annotation} \
            {output.json} \
            2> {log}
        """


rule tair_to_go:
    """
    Map TAIR Arabidopsis gene IDs to GO terms.

    Uses TAIR gene association file (GAF format).
    """
    input:
        annotation = ANNOTATION_INFO,
        tair_assoc = TAIR_ASSOC
    output:
        json = GO_MAPPING_DIR / "tair_to_go.json"
    log:
        GO_LOG_DIR / "tair_to_go.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/go_analysis/tair_to_go.py \
            {input.annotation} \
            {input.tair_assoc} \
            {output.json} \
            2> {log}
        """


rule panther_to_go:
    """
    Parse PANTHER database for GO term mappings.
    """
    input:
        panther_db = PANTHER_HMM
    output:
        json = GO_MAPPING_DIR / "panther_to_go.json"
    log:
        GO_LOG_DIR / "panther_to_go.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/go_analysis/panther_to_go.py \
            {input.panther_db} \
            {output.json} \
            2> {log}
        """


rule enhance_go_coverage:
    """
    Merge all GO mappings into comprehensive gene-to-GO dictionary.

    Combines:
    - Direct GO annotations from annotation file
    - Pfam-derived GO terms
    - KO-derived GO terms
    - TAIR-derived GO terms
    - PANTHER-derived GO terms
    """
    input:
        annotation = ANNOTATION_INFO,
        pfam_json = GO_MAPPING_DIR / "pfam_to_go.json",
        ko_json = GO_MAPPING_DIR / "ko_to_go.json",
        tair_json = GO_MAPPING_DIR / "tair_to_go.json",
        panther_json = GO_MAPPING_DIR / "panther_to_go.json"
    output:
        json = GO_MAPPING_DIR / "gene2go.json",
        stats = GO_MAPPING_DIR / "go_coverage_stats.txt"
    log:
        GO_LOG_DIR / "enhance_go_coverage.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/go_analysis/enhance_go_coverage.py \
            {input.annotation} \
            {input.pfam_json} \
            {input.ko_json} \
            {input.tair_json} \
            {input.panther_json} \
            {output.json} \
            2> {log}

        # Generate coverage statistics
        python -c "
import json
with open('{output.json}') as f:
    data = json.load(f)
n_genes = len(data)
n_with_go = sum(1 for g in data.values() if g)
total_go = sum(len(g) for g in data.values())
with open('{output.stats}', 'w') as f:
    f.write('GO Coverage Statistics\\n')
    f.write('=' * 40 + '\\n')
    f.write(f'Total genes: {{n_genes}}\\n')
    f.write(f'Genes with GO terms: {{n_with_go}} ({{100*n_with_go/n_genes:.1f}}%)\\n')
    f.write(f'Total GO annotations: {{total_go}}\\n')
    f.write(f'Average GO terms per gene: {{total_go/n_genes:.1f}}\\n')
"
        """


# ============================================================
# PHASE 1: GENE CLASSIFICATION BY INTRONER STATUS
# ============================================================

rule introner_phase1_analysis:
    """
    Classify genes by introner presence patterns.

    For each gene, determines:
    - group1_fixed: Introner fixed in all Group1 samples
    - group1_polymorphic: Introner polymorphic within Group1
    - group2_fixed: Introner fixed in all Group2 samples
    - group2_polymorphic: Introner polymorphic within Group2
    """
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.tsv",
        gtf_dir = ANNOTATIONS_DIR,
        go_json = GO_MAPPING_DIR / "gene2go.json"
    output:
        results = GO_RESULTS_DIR / "phase1_gene_classification.tsv",
        summary = GO_RESULTS_DIR / "phase1_summary.txt"
    log:
        GO_LOG_DIR / "phase1_analysis.log"
    shell:
        """
        mkdir -p {GO_RESULTS_DIR}

        python {PROJECT_ROOT}/scripts/go_analysis/introner_phase1_analysis.py \
            {input.genotype_matrix} \
            {input.gtf_dir} \
            {input.go_json} \
            {output.results} \
            2> {log}

        # Generate summary
        python -c "
import pandas as pd
df = pd.read_csv('{output.results}', sep='\\t')
with open('{output.summary}', 'w') as f:
    f.write('Phase 1 Gene Classification Summary\\n')
    f.write('=' * 50 + '\\n\\n')
    f.write(f'Total genes analyzed: {{len(df)}}\\n')
    f.write(f'Genes with Group1 fixed introners: {{df[\"group1_fixed\"].sum()}}\\n')
    f.write(f'Genes with Group1 polymorphic introners: {{df[\"group1_polymorphic\"].sum()}}\\n')
    f.write(f'Genes with Group2 fixed introners: {{df[\"group2_fixed\"].sum()}}\\n')
    f.write(f'Genes with Group2 polymorphic introners: {{df[\"group2_polymorphic\"].sum()}}\\n')
    f.write(f'Genes with GO annotations: {{(df[\"go_terms_str\"].str.len() > 0).sum()}}\\n')
"
        """


rule plot_phase1_results:
    """
    Visualize Phase 1 gene classification results.
    """
    input:
        results = GO_RESULTS_DIR / "phase1_gene_classification.tsv",
        go_json = GO_MAPPING_DIR / "gene2go.json"
    output:
        overview = GO_PLOTS_DIR / "phase1_overview.pdf",
        venn = GO_PLOTS_DIR / "phase1_venn.pdf"
    log:
        GO_LOG_DIR / "plot_phase1.log"
    shell:
        """
        mkdir -p {GO_PLOTS_DIR}

        python {PROJECT_ROOT}/scripts/go_analysis/plot_phase1_results.py \
            {input.results} \
            {GO_PLOTS_DIR} \
            --go_json {input.go_json} \
            2> {log}
        """


# ============================================================
# PHASE 2: GO ENRICHMENT ANALYSIS
# ============================================================

rule introner_phase2_enrichment:
    """
    Perform GO enrichment analysis for introner categories.

    Tests for functional enrichment in:
    - Group1 fixed introner genes
    - Group1 polymorphic introner genes
    - Group2 fixed introner genes
    - Group2 polymorphic introner genes

    Uses Fisher's exact test with FDR correction.
    """
    input:
        phase1 = GO_RESULTS_DIR / "phase1_gene_classification.tsv",
        go_json = GO_MAPPING_DIR / "gene2go.json"
    output:
        results = GO_RESULTS_DIR / "phase2_enrichment_results.tsv",
        significant = GO_RESULTS_DIR / "phase2_significant_terms.tsv"
    params:
        fdr_threshold = config["params"]["go"]["fdr_threshold"],
        min_genes = config["params"]["go"]["min_genes"]
    log:
        GO_LOG_DIR / "phase2_enrichment.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/go_analysis/introner_phase2_enrichment.py \
            {input.phase1} \
            {output.results} \
            --go_json {input.go_json} \
            --fdr_threshold {params.fdr_threshold} \
            --min_genes {params.min_genes} \
            2> {log}

        # Extract significant terms
        python -c "
import pandas as pd
df = pd.read_csv('{output.results}', sep='\\t')
sig = df[df['significant'] == True] if 'significant' in df.columns else df[df['fdr_pvalue'] < {params.fdr_threshold}]
sig.to_csv('{output.significant}', sep='\\t', index=False)
"
        """


rule introner_phase2_family_enrichment:
    """
    Perform family-stratified GO enrichment analysis.

    Analyzes GO enrichment separately for each introner family type.
    """
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.tsv",
        phase1 = GO_RESULTS_DIR / "phase1_gene_classification.tsv",
        go_json = GO_MAPPING_DIR / "gene2go.json"
    output:
        results = GO_RESULTS_DIR / "phase2_family_enrichment_results.tsv"
    log:
        GO_LOG_DIR / "phase2_family_enrichment.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/go_analysis/introner_phase2_family_enrichment.py \
            {input.genotype_matrix} \
            {input.phase1} \
            {output.results} \
            --go_json {input.go_json} \
            2> {log}
        """


rule plot_phase2_results:
    """
    Visualize Phase 2 GO enrichment results.

    Creates:
    - Enrichment heatmap (GO categories x introner types)
    - Volcano plots (effect size vs significance)
    - Top enriched categories bar charts
    """
    input:
        results = GO_RESULTS_DIR / "phase2_enrichment_results.tsv"
    output:
        heatmap = GO_PLOTS_DIR / "phase2_enrichment_heatmap.pdf",
        volcano = GO_PLOTS_DIR / "phase2_volcano_plots.pdf",
        barchart = GO_PLOTS_DIR / "phase2_top_enriched.pdf"
    log:
        GO_LOG_DIR / "plot_phase2.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/go_analysis/plot_phase2_results.py \
            {input.results} \
            {GO_PLOTS_DIR} \
            2> {log}
        """


# ============================================================
# STANDARD GO ENRICHMENT (GOATOOLS-BASED)
# ============================================================

rule run_goatools_enrichment:
    """
    Run standard GO enrichment using goatools library.

    Provides additional statistical methods:
    - Bonferroni correction
    - Sidak correction
    - Holm correction
    - FDR (Benjamini-Hochberg)
    """
    input:
        phase1 = GO_RESULTS_DIR / "phase1_gene_classification.tsv",
        go_obo = PROJECT_ROOT / "resources" / "go.obo",
        go_json = GO_MAPPING_DIR / "gene2go.json"
    output:
        enrichment = GO_RESULTS_DIR / "goatools_enrichment_results.tsv"
    log:
        GO_LOG_DIR / "goatools_enrichment.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/go_analysis/go_enrichment.py \
            --go_obo {input.go_obo} \
            --go_json {input.go_json} \
            --genes {input.phase1} \
            --output {output.enrichment} \
            2> {log}
        """


# ============================================================
# COPY FIGURES TO MAIN FIGURES DIRECTORY
# ============================================================

rule copy_go_figures:
    """
    Copy key GO analysis figures to main figures directory.
    """
    input:
        heatmap = GO_PLOTS_DIR / "phase2_enrichment_heatmap.pdf",
        overview = GO_PLOTS_DIR / "phase1_overview.pdf"
    output:
        heatmap = FIGURES_DIR / "go_enrichment_heatmap.pdf",
        overview = FIGURES_DIR / "go_phase1_overview.pdf"
    shell:
        """
        cp {input.heatmap} {output.heatmap}
        cp {input.overview} {output.overview}
        """


# ============================================================
# TARGET RULES
# ============================================================

rule go_enrichment_complete:
    """
    Target: Complete GO enrichment analysis pipeline.
    """
    input:
        # GO mappings
        GO_MAPPING_DIR / "gene2go.json",
        GO_MAPPING_DIR / "go_coverage_stats.txt",
        # Phase 1
        GO_RESULTS_DIR / "phase1_gene_classification.tsv",
        GO_PLOTS_DIR / "phase1_overview.pdf",
        # Phase 2
        GO_RESULTS_DIR / "phase2_enrichment_results.tsv",
        GO_RESULTS_DIR / "phase2_family_enrichment_results.tsv",
        GO_PLOTS_DIR / "phase2_enrichment_heatmap.pdf",
        # Main figures
        FIGURES_DIR / "go_enrichment_heatmap.pdf"


rule go_mapping_only:
    """
    Target: Generate GO mappings only.
    """
    input:
        GO_MAPPING_DIR / "gene2go.json",
        GO_MAPPING_DIR / "go_coverage_stats.txt"


rule go_phase1_only:
    """
    Target: Phase 1 gene classification only.
    """
    input:
        GO_RESULTS_DIR / "phase1_gene_classification.tsv",
        GO_PLOTS_DIR / "phase1_overview.pdf"


rule go_phase2_only:
    """
    Target: Phase 2 enrichment analysis only.
    """
    input:
        GO_RESULTS_DIR / "phase2_enrichment_results.tsv",
        GO_PLOTS_DIR / "phase2_enrichment_heatmap.pdf"
