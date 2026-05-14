# ============================================================
# 40_go_enrichment.smk - GO Term Enrichment Analysis
# ============================================================
#
# Performs Gene Ontology enrichment analysis for biologically defined
# introner groups in genotype_matrix.final.tsv:
# 1. Generate GO mappings from multiple databases (Pfam, KO, TAIR, PANTHER)
# 2. Enhance GO coverage by merging all sources
# 3. Classify final-matrix introner loci as all introners, ancestral,
#    independent insertion, polymorphic within Group 1, or consistently
#    present within Group 1/2
# 4. Run full-GO and GO Slim enrichment and generate final figures
#
# Adapted from:
# - /scratch1/chris/mpusilla_go_analysis/scripts/
#
# ============================================================

from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Output directories
GO_DIR = RESULTS / "go_enrichment"
GO_MAPPING_DIR = GO_DIR / "mappings"
GO_RESULTS_DIR = GO_DIR / "results"
GO_LOG_DIR = GO_DIR / "logs"

# External database paths
INTERPRO2GO = config["paths"]["external_db"]["interpro2go"]
INTERPRO_TO_PFAM = config["paths"]["external_db"].get(
    "interpro_to_pfam",
    str(Path(INTERPRO2GO).parent / "interpro_to_pfam.json")
)
TAIR_ASSOC = config["paths"]["external_db"]["gene_association_tair"]
PANTHER_HMM = config["paths"]["external_db"]["panther_hmm"]

# Annotation file
ANNOTATION_INFO = config["paths"]["references"]["annotation_info"]
GO_API_TIMEOUT = config["params"]["go"].get("api_timeout", 30)
GO_API_MAX_RETRIES = config["params"]["go"].get("api_max_retries", 3)


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
        interpro_to_pfam = INTERPRO_TO_PFAM,
        annotation = ANNOTATION_INFO
    output:
        json = GO_MAPPING_DIR / "pfam_to_go.json"
    params:
        timeout = GO_API_TIMEOUT,
        max_retries = GO_API_MAX_RETRIES
    log:
        GO_LOG_DIR / "pfam_to_go.log"
    conda:
        "../envs/go_enrichment.yaml"
    shell:
        """
        mkdir -p {GO_MAPPING_DIR}
        mkdir -p {GO_LOG_DIR}

        python {PROJECT_ROOT}/scripts/go_analysis/pfam_to_go.py \
            {input.interpro2go} \
            {input.annotation} \
            {output.json} \
            --interpro-to-pfam {input.interpro_to_pfam} \
            --timeout {params.timeout} \
            --max-retries {params.max_retries} \
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
    params:
        batch_size = config["params"]["go"].get("kegg_batch_size", 10),
        timeout = GO_API_TIMEOUT,
        max_retries = GO_API_MAX_RETRIES
    log:
        GO_LOG_DIR / "ko_to_go.log"
    conda:
        "../envs/go_enrichment.yaml"
    shell:
        """
        mkdir -p {GO_MAPPING_DIR}
        mkdir -p {GO_LOG_DIR}

        python {PROJECT_ROOT}/scripts/go_analysis/ko_to_go.py \
            {input.annotation} \
            {output.json} \
            --batch-size {params.batch_size} \
            --timeout {params.timeout} \
            --max-retries {params.max_retries} \
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
    conda:
        "../envs/go_enrichment.yaml"
    shell:
        """
        mkdir -p {GO_MAPPING_DIR}
        mkdir -p {GO_LOG_DIR}

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
    conda:
        "../envs/go_enrichment.yaml"
    shell:
        """
        mkdir -p {GO_MAPPING_DIR}
        mkdir -p {GO_LOG_DIR}

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
    conda:
        "../envs/go_enrichment.yaml"
    shell:
        """
        mkdir -p {GO_MAPPING_DIR}
        mkdir -p {GO_LOG_DIR}

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
    pct_with_go = 100*n_with_go/n_genes if n_genes else 0.0
    avg_go = total_go/n_genes if n_genes else 0.0
    f.write(f'Genes with GO terms: {{n_with_go}} ({{pct_with_go:.1f}}%)\\n')
    f.write(f'Total GO annotations: {{total_go}}\\n')
    f.write(f'Average GO terms per gene: {{avg_go:.1f}}\\n')
"
        """


# ============================================================
# FINAL-MATRIX INTRONER GROUP GO ENRICHMENT
# ============================================================

rule introner_group_go_enrichment:
    """
    GO enrichment for biologically defined introner groups in genotype_matrix.final.tsv.

    Groups:
    - ancestral: cross_group_status in ancestral/likely_ancestral/ancestral_low_identity
    - independent_insertion: cross_group_status in independent/likely_independent
    - polymorphic_within_group1: presence=1 and presence=2 both observed within Group 1
    - consistent_group1: presence=1 for every Group 1 sample
    - consistent_group2: presence=1 for every Group 2 sample
    - all_introners: presence=1 in at least one Group 1 or Group 2 sample

    Missing calls (presence=3) are ignored for polymorphic and cross-group present
    evidence. Consistent groups require every configured sample in that group to
    be present. Background is all GO-annotated reference genes by default,
    which supports the all_introners contrast against genes without introners.
    """
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        go_json = GO_MAPPING_DIR / "gene2go.json",
        go_obo = PROJECT_ROOT / "resources" / "go.obo",
        script = PROJECT_ROOT / "scripts" / "go_analysis" / "introner_group_go_enrichment.py"
    output:
        loci = GO_RESULTS_DIR / "introner_group_locus_classification.tsv",
        gene_sets = GO_RESULTS_DIR / "introner_group_gene_sets.tsv",
        enrichment = GO_RESULTS_DIR / "introner_group_go_enrichment.tsv",
        significant = GO_RESULTS_DIR / "introner_group_go_enrichment.significant.tsv",
        summary = GO_RESULTS_DIR / "introner_group_go_enrichment.summary.txt"
    params:
        group1_samples = ",".join(GROUP1_SAMPLES),
        group2_samples = ",".join(GROUP2_SAMPLES),
        fdr_threshold = config["params"]["go"]["fdr_threshold"],
        min_genes = config["params"]["go"]["min_genes"]
    log:
        GO_LOG_DIR / "introner_group_go_enrichment.log"
    conda:
        "../envs/go_enrichment.yaml"
    shell:
        """
        mkdir -p {GO_RESULTS_DIR}
        mkdir -p {GO_LOG_DIR}

        python {input.script} \
            --genotype-matrix {input.genotype_matrix} \
            --go-json {input.go_json} \
            --go-obo {input.go_obo} \
            --group1-samples {params.group1_samples} \
            --group2-samples {params.group2_samples} \
            --locus-output {output.loci} \
            --gene-set-output {output.gene_sets} \
            --enrichment-output {output.enrichment} \
            --significant-output {output.significant} \
            --summary-output {output.summary} \
            --fdr-threshold {params.fdr_threshold} \
            --min-genes {params.min_genes} \
            --background-scope all_annotated \
            2> {log}
        """


rule introner_group_go_slim_enrichment:
    """
    GO Slim enrichment for the same final-matrix introner groups.

    Uses broad propagated GO ancestor categories to reduce the multiple-testing
    burden and test general enrichment/depletion patterns.
    """
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        go_json = GO_MAPPING_DIR / "gene2go.json",
        go_obo = PROJECT_ROOT / "resources" / "go.obo",
        script = PROJECT_ROOT / "scripts" / "go_analysis" / "introner_group_go_slim_enrichment.py",
        classifier = PROJECT_ROOT / "scripts" / "go_analysis" / "introner_group_go_enrichment.py"
    output:
        enrichment = GO_RESULTS_DIR / "introner_group_go_slim_enrichment.tsv",
        significant = GO_RESULTS_DIR / "introner_group_go_slim_enrichment.significant.tsv",
        summary = GO_RESULTS_DIR / "introner_group_go_slim_enrichment.summary.txt"
    params:
        group1_samples = ",".join(GROUP1_SAMPLES),
        group2_samples = ",".join(GROUP2_SAMPLES),
        fdr_threshold = config["params"]["go"]["fdr_threshold"],
        min_genes = config["params"]["go"]["min_genes"]
    log:
        GO_LOG_DIR / "introner_group_go_slim_enrichment.log"
    conda:
        "../envs/go_enrichment.yaml"
    shell:
        """
        mkdir -p {GO_RESULTS_DIR}
        mkdir -p {GO_LOG_DIR}

        python {input.script} \
            --genotype-matrix {input.genotype_matrix} \
            --go-json {input.go_json} \
            --go-obo {input.go_obo} \
            --group1-samples {params.group1_samples} \
            --group2-samples {params.group2_samples} \
            --enrichment-output {output.enrichment} \
            --significant-output {output.significant} \
            --summary-output {output.summary} \
            --fdr-threshold {params.fdr_threshold} \
            --min-genes {params.min_genes} \
            --background-scope all_annotated \
            2> {log}
        """


rule plot_introner_group_go_enrichment:
    """
    Plot GO enrichment results generated from genotype_matrix.final.tsv.
    """
    input:
        enrichment = GO_RESULTS_DIR / "introner_group_go_enrichment.tsv",
        script = PROJECT_ROOT / "scripts" / "go_analysis" / "plot_introner_group_go_enrichment.py"
    output:
        heatmap = FIGURES_DIR / "go_enrichment_heatmap.pdf",
        heatmap_png = FIGURES_DIR / "go_enrichment_heatmap.png",
        top_terms = FIGURES_DIR / "go_enrichment_top_terms.pdf",
        top_terms_png = FIGURES_DIR / "go_enrichment_top_terms.png"
    log:
        GO_LOG_DIR / "plot_introner_group_go_enrichment.log"
    conda:
        "../envs/go_enrichment.yaml"
    shell:
        """
        mkdir -p {FIGURES_DIR}
        mkdir -p {GO_LOG_DIR}
        mkdir -p {GO_LOG_DIR}/matplotlib

        MPLCONFIGDIR={GO_LOG_DIR}/matplotlib \
            python {input.script} \
            {input.enrichment} \
            --heatmap-output {output.heatmap} \
            --top-terms-output {output.top_terms} \
            2> {log}
        """


# ============================================================
# TARGET RULES
# ============================================================

rule go_enrichment_complete:
    """
    Target: Complete current GO enrichment analysis pipeline.

    Uses the final genotype matrix and the biologically defined introner groups:
    all_introners, ancestral, independent_insertion, polymorphic_within_group1,
    consistent_group1, and consistent_group2.
    """
    input:
        # GO mappings
        GO_MAPPING_DIR / "gene2go.json",
        GO_MAPPING_DIR / "go_coverage_stats.txt",
        # Final-matrix introner group enrichment
        GO_RESULTS_DIR / "introner_group_locus_classification.tsv",
        GO_RESULTS_DIR / "introner_group_gene_sets.tsv",
        GO_RESULTS_DIR / "introner_group_go_enrichment.tsv",
        GO_RESULTS_DIR / "introner_group_go_enrichment.significant.tsv",
        GO_RESULTS_DIR / "introner_group_go_enrichment.summary.txt",
        GO_RESULTS_DIR / "introner_group_go_slim_enrichment.tsv",
        GO_RESULTS_DIR / "introner_group_go_slim_enrichment.significant.tsv",
        GO_RESULTS_DIR / "introner_group_go_slim_enrichment.summary.txt",
        FIGURES_DIR / "go_enrichment_heatmap.pdf",
        FIGURES_DIR / "go_enrichment_top_terms.pdf"


rule go_mapping_only:
    """
    Target: Generate GO mappings only.
    """
    input:
        GO_MAPPING_DIR / "gene2go.json",
        GO_MAPPING_DIR / "go_coverage_stats.txt"


rule introner_group_go_enrichment_complete:
    """
    Target: Final-matrix GO enrichment for all introners, ancestral,
    independent insertion, Group 1 polymorphic, and consistently present
    introner groups.
    """
    input:
        GO_RESULTS_DIR / "introner_group_locus_classification.tsv",
        GO_RESULTS_DIR / "introner_group_gene_sets.tsv",
        GO_RESULTS_DIR / "introner_group_go_enrichment.tsv",
        GO_RESULTS_DIR / "introner_group_go_enrichment.significant.tsv",
        GO_RESULTS_DIR / "introner_group_go_enrichment.summary.txt",
        GO_RESULTS_DIR / "introner_group_go_slim_enrichment.tsv",
        GO_RESULTS_DIR / "introner_group_go_slim_enrichment.significant.tsv",
        GO_RESULTS_DIR / "introner_group_go_slim_enrichment.summary.txt",
        FIGURES_DIR / "go_enrichment_heatmap.pdf",
        FIGURES_DIR / "go_enrichment_top_terms.pdf"
