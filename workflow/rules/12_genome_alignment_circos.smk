# ============================================================================
# 12_genome_alignment_circos.smk
# Whole-genome MUMmer alignment and Circos synteny visualization
# Compares CCMP1545 (intronerful) vs RCC1749 (intronerless) assemblies
#
# Source: /scratch1/chris/introner-group-comparison-analysis/rules/mummer_circos.smk
# ============================================================================

import os

# Strain identifiers for the comparison
G1 = config["samples"]["reference"]["group1"]   # CCMP1545
G2 = config["samples"]["reference"]["group2"]   # RCC1749

# Output subdirectories
MUMMER_DIR = GENOME_ALIGNMENT_DIR / "mummer"
CIRCOS_DIR = GENOME_ALIGNMENT_DIR / "circos"
GA_FIGURES_DIR = FIGURES_DIR / "genome_alignment"

# Script directory
GA_SCRIPTS = str(Path(workflow.basedir).parent / "scripts" / "genome_alignment")

# Config params
MUMMER_PARAMS = config.get("params", {}).get("mummer", {})
CIRCOS_PARAMS = config.get("params", {}).get("circos", {})

# ============================================================================
# MUMMER ALIGNMENT RULES
# ============================================================================

rule nucmer_alignment:
    """Run nucmer whole-genome alignment and filter for 1-to-1 mappings."""
    input:
        g1_genome = ASSEMBLIES_DIR / f"{G1}.vg_paths.fa",
        g2_genome = ASSEMBLIES_DIR / f"{G2}.vg_paths.fa"
    output:
        delta = MUMMER_DIR / f"{G1}_vs_{G2}.1delta"
    params:
        prefix = str(MUMMER_DIR / f"{G1}_vs_{G2}")
    shell:
        """
        mkdir -p {MUMMER_DIR} && \
        nucmer {input.g1_genome} {input.g2_genome} -p {params.prefix} && \
        delta-filter -1 {params.prefix}.delta > {output.delta}
        """

rule mummerplot_dotplot:
    """Generate dot plot visualization from MUMmer alignment."""
    input:
        delta = MUMMER_DIR / f"{G1}_vs_{G2}.1delta"
    output:
        png = GA_FIGURES_DIR / f"{G1}_vs_{G2}.png"
    params:
        prefix = str(GA_FIGURES_DIR / f"{G1}_vs_{G2}")
    shell:
        """
        mkdir -p {GA_FIGURES_DIR} && \
        mummerplot --png -p {params.prefix} --large {input.delta}
        """

# ============================================================================
# CIRCOS PREPARATION RULES
# ============================================================================

rule generate_karyotype:
    """Generate circos karyotype file from assembly indices."""
    input:
        g1_idx = ASSEMBLIES_DIR / f"{G1}.vg_paths.fa.fai",
        g2_idx = ASSEMBLIES_DIR / f"{G2}.vg_paths.fa.fai"
    output:
        kar = CIRCOS_DIR / f"{G1}_vs_{G2}.kar"
    shell:
        """
        mkdir -p {CIRCOS_DIR} && \
        python3 {GA_SCRIPTS}/generate_karyotype.py \
            --strain1-fai {input.g1_idx} \
            --strain2-fai {input.g2_idx} \
            --strain1-name {G1} \
            --strain2-name {G2} \
            --output {output.kar}
        """

rule delta_to_coords:
    """Convert filtered delta file to show-coords format for link generation."""
    input:
        delta = MUMMER_DIR / f"{G1}_vs_{G2}.1delta"
    output:
        coords = MUMMER_DIR / f"{G1}_vs_{G2}.coords"
    shell:
        """
        show-coords -c -l -T {input.delta} > {output.coords}
        """

rule generate_links:
    """Extract synteny blocks from MUMmer coords file."""
    input:
        coords = MUMMER_DIR / f"{G1}_vs_{G2}.coords"
    output:
        links = CIRCOS_DIR / f"{G1}_vs_{G2}.links.tsv"
    params:
        min_length = MUMMER_PARAMS.get("min_link_length", 1000),
        min_identity = MUMMER_PARAMS.get("min_identity", 80.0)
    shell:
        """
        python3 {GA_SCRIPTS}/generate_links.py \
            --coords {input.coords} \
            --strain1-name {G1} \
            --strain2-name {G2} \
            --output {output.links} \
            --min-length {params.min_length} \
            --min-identity {params.min_identity}
        """

rule analyze_synteny:
    """Reorder karyotype based on synteny relationships."""
    input:
        links = CIRCOS_DIR / f"{G1}_vs_{G2}.links.tsv",
        karyotype = CIRCOS_DIR / f"{G1}_vs_{G2}.kar"
    output:
        kar = CIRCOS_DIR / f"{G1}_vs_{G2}_reordered.kar"
    shell:
        """
        python3 {GA_SCRIPTS}/analyze_synteny.py \
            --links {input.links} \
            --karyotype {input.karyotype} \
            --output-karyotype {output.kar}
        """

rule flip_chromosomes:
    """Flip RCC1749 contig orientations to match CCMP1545."""
    input:
        kar = CIRCOS_DIR / f"{G1}_vs_{G2}_reordered.kar"
    output:
        kar = CIRCOS_DIR / f"{G1}_vs_{G2}_reordered_flipped.kar"
    shell:
        """
        python3 {GA_SCRIPTS}/flip_rcc1749_order.py \
            --input {input.kar} \
            --output {output.kar}
        """

rule filter_reorder_karyotype:
    """Remove small/uninformative contigs and adjust contig ordering."""
    input:
        kar = CIRCOS_DIR / f"{G1}_vs_{G2}_reordered_flipped.kar"
    output:
        kar = CIRCOS_DIR / f"{G1}_vs_{G2}_reordered_flipped_filtered.kar"
    params:
        remove_contigs = CIRCOS_PARAMS.get("remove_contigs", "2,1,33,30,32,40,43,24,42"),
        move_contig = CIRCOS_PARAMS.get("move_contig", 20),
        insert_after = CIRCOS_PARAMS.get("insert_after", 25),
        insert_before = CIRCOS_PARAMS.get("insert_before", 12),
        move_contig_color = CIRCOS_PARAMS.get("move_contig_color", "")
    shell:
        """
        python3 {GA_SCRIPTS}/filter_reorder_karyotype.py \
            --input {input.kar} \
            --output {output.kar} \
            --remove-contigs {params.remove_contigs} \
            --move-contig {params.move_contig} \
            --insert-after {params.insert_after} \
            --insert-before {params.insert_before} \
            --move-contig-color {params.move_contig_color}
        """

rule generate_ribbons:
    """Generate colored synteny ribbons for circos visualization."""
    input:
        karyotype = CIRCOS_DIR / f"{G1}_vs_{G2}_reordered_flipped_filtered.kar",
        links = CIRCOS_DIR / f"{G1}_vs_{G2}.links.tsv"
    output:
        ribbons = CIRCOS_DIR / f"{G1}_vs_{G2}_ribbons_with_merging.tsv"
    params:
        gap_threshold = CIRCOS_PARAMS.get("gap_threshold", 50000),
        merge_threshold = CIRCOS_PARAMS.get("merge_threshold", 750),
        min_link_length = MUMMER_PARAMS.get("min_link_length", 1000)
    shell:
        """
        python3 {GA_SCRIPTS}/generate_ribbon_links_with_colors.py \
            --karyotype {input.karyotype} \
            --input {input.links} \
            --gap-threshold {params.gap_threshold} \
            --min-link-length {params.min_link_length} \
            --merge-threshold {params.merge_threshold} \
            --output {output.ribbons}
        """

# ============================================================================
# CIRCOS HISTOGRAM RULES
# ============================================================================

# Non-introner intron matrix is produced by rules/13_non_introner_introns.smk.
NON_INTRONER_INTRONS_DIR_FOR_CIRCOS = EVOLUTION_DIR / "non_introner_introns"

rule generate_circos_histograms:
    """Generate Circos histogram tracks from current genotype matrices."""
    input:
        karyotype = CIRCOS_DIR / f"{G1}_vs_{G2}_reordered_flipped_filtered.kar",
        introner_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        non_introner_matrix = NON_INTRONER_INTRONS_DIR_FOR_CIRCOS / "intron_genotype_matrix.tsv",
        g1_gtf = get_gtf(G1),
        g2_gtf = get_gtf(G2)
    output:
        g1_introns = CIRCOS_DIR / f"{G1}.intron_histogram.txt",
        g2_introns = CIRCOS_DIR / f"{G2}.intron_histogram.txt",
        g1_introners = CIRCOS_DIR / f"{G1}.introner_histogram.txt",
        g2_introners = CIRCOS_DIR / f"{G2}.introner_histogram.txt"
    params:
        bin_size = CIRCOS_PARAMS.get("histogram_bin_size", 50000),
        output_dir = CIRCOS_DIR
    shell:
        """
        python3 {GA_SCRIPTS}/generate_circos_histograms.py \
            --introner-matrix {input.introner_matrix} \
            --non-introner-matrix {input.non_introner_matrix} \
            --karyotype {input.karyotype} \
            --strain1-name {G1} \
            --strain2-name {G2} \
            --non-introner-reference {REFERENCE} \
            --strain1-gtf {input.g1_gtf} \
            --strain2-gtf {input.g2_gtf} \
            --bin-size {params.bin_size} \
            --output-dir {params.output_dir}
        """

# ============================================================================
# CIRCOS RENDERING RULES
# ============================================================================

# Static circos config lives in data/circos/.
CIRCOS_DATA_DIR = PROJECT_ROOT / "data" / "circos"

rule render_circos:
    """Render final circos synteny plot using generated data tracks."""
    input:
        karyotype = CIRCOS_DIR / f"{G1}_vs_{G2}_reordered_flipped_filtered.kar",
        ribbons = CIRCOS_DIR / f"{G1}_vs_{G2}_ribbons_with_merging.tsv",
        conf = CIRCOS_DATA_DIR / "circos.conf",
        histograms = rules.generate_circos_histograms.output
    output:
        png = CIRCOS_DIR / "circos.png",
        svg = CIRCOS_DIR / "circos.svg"
    params:
        circos_env = config.get("params", {}).get("circos", {}).get("conda_env", "circos_env"),
        data_dir = CIRCOS_DATA_DIR
    shell:
        """
        cp {params.data_dir}/circos.conf {CIRCOS_DIR}/ && \
        cd {CIRCOS_DIR} && conda run -n {params.circos_env} circos -conf circos.conf
        """
