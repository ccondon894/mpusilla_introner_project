# Common exonic backgrounds, matched class comparisons and SNP enrichment (38–39).

rule build_common_recombination_background:
    input:
        script = PROJECT_ROOT / "scripts/popgen/recombination_analysis/build_common_recombination_background.py",
        all = SNP_DIR / "selection/recombination" / "gene_exonic_introners_all_5kb_updated_windows.tsv",
        poly = SNP_DIR / "selection/recombination" / "gene_exonic_polymorphic_introners_5kb_updated_windows.tsv",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        canonical = EVOLUTION_DIR / "non_introner_introns/intron_genotype_matrix.tsv",
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        maps = sorted(Path(config["paths"]["pyrho"]["ccmp1545"]).glob("*.pyrho.out"))
    output:
        file0 = SNP_DIR / "selection/recombination" / "common_background" / "group1.background.tsv",
        file1 = SNP_DIR / "selection/recombination" / "common_background" / "group1.excluded_background.tsv",
        file2 = SNP_DIR / "selection/recombination" / "common_background" / "group1.fixed_focals.tsv",
        file3 = SNP_DIR / "selection/recombination" / "common_background" / "group1.canonical_focals.tsv",
        file4 = SNP_DIR / "selection/recombination" / "common_background" / "group1.comparison_windows.tsv",
        file5 = SNP_DIR / "selection/recombination" / "common_background" / "group1.audit.tsv"
    params:
        outdir = SNP_DIR / "selection/recombination" / "common_background",
        pyrho = config["paths"]["pyrho"]["ccmp1545"]
    threads: 1
    log:
        RESULTS / "logs/promoted_analyses/build_common_recombination_background.log"
    conda: "../envs/analysis_tools.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --all-introner-windows {input.all:q} --polymorphic-introner-windows {input.poly:q} --introner-genotype-matrix {input.matrix:q} --canonical-intron-matrix {input.canonical:q} --gtf {input.gtf:q} --pyrho-dir {params.pyrho:q} --output-background {output.file0:q} --output-excluded-background {output.file1:q} --output-fixed-introner-focals {output.file2:q} --output-canonical-focals {output.file3:q} --output-comparison-windows {output.file4:q} --output-audit {output.file5:q} > {log:q} 2>&1
        """

rule build_population2_recombination_background:
    input:
        script = PROJECT_ROOT / "scripts/popgen/recombination_analysis/build_population2_background.py",
        windows = SNP_DIR / "selection/recombination" / "gene_exonic_group2_introners_5kb_updated_windows.tsv"
    output:
        file0 = SNP_DIR / "selection/recombination" / "common_background" / "group2.background.tsv",
        file1 = SNP_DIR / "selection/recombination" / "common_background" / "group2.excluded_background.tsv",
        file2 = SNP_DIR / "selection/recombination" / "common_background" / "group2.focals.tsv",
        file3 = SNP_DIR / "selection/recombination" / "common_background" / "group2.comparison_windows.tsv",
        file4 = SNP_DIR / "selection/recombination" / "common_background" / "group2.audit.tsv"
    params:
        outdir = SNP_DIR / "selection/recombination" / "common_background"
    threads: 1
    log:
        RESULTS / "logs/promoted_analyses/build_population2_recombination_background.log"
    conda: "../envs/analysis_tools.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --group2-windows {input.windows:q} --output-background {output.file0:q} --output-excluded-background {output.file1:q} --output-focals {output.file2:q} --output-comparison-windows {output.file3:q} --output-audit {output.file4:q} > {log:q} 2>&1
        """

rule recombination_common_background_tests:
    input:
        script = PROJECT_ROOT / "scripts/popgen/recombination_analysis/run_mann_whitney_comparisons.py",
        group1 = SNP_DIR / "selection/recombination" / "common_background" / "group1.comparison_windows.tsv",
        group2 = SNP_DIR / "selection/recombination" / "common_background" / "group2.comparison_windows.tsv"
    output:
        file0 = SNP_DIR / "selection/recombination" / "common_background" / "mann_whitney_results.tsv"
    params:
        outdir = SNP_DIR / "selection/recombination" / "common_background"
    threads: 1
    log:
        RESULTS / "logs/promoted_analyses/recombination_common_background_tests.log"
    conda: "../envs/analysis_tools.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --group1-comparison-windows {input.group1:q} --group2-windows {input.group2:q} --output {output.file0:q} > {log:q} 2>&1
        """

rule snp_recombination_enrichment:
    input:
        script = PROJECT_ROOT / "scripts/popgen/recombination_analysis/analyze_snp_recombination_enrichment.py",
        vcf = config["paths"]["fourfold_vcf"],
        maps = sorted(Path(config["paths"]["pyrho"]["ccmp1545"]).glob("*.pyrho.out")),
        colors = PROJECT_ROOT / "master_figure_color_guide.tsv",
        color_script = PROJECT_ROOT / "scripts/figure_color_guide.py"
    output:
        report = directory(SNP_DIR / "selection/snp_recombination_enrichment")
    params:
        outdir = SNP_DIR / "selection/snp_recombination_enrichment",
        pyrho = config["paths"]["pyrho"]["ccmp1545"]
    threads: 1
    log:
        RESULTS / "logs/promoted_analyses/snp_recombination_enrichment.log"
    conda: "../envs/analysis_tools.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --vcf {input.vcf:q} --pyrho-dir {params.pyrho:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """
