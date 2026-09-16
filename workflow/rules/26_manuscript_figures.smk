# Manuscript and presentation figure assemblies (54–57).

rule introner_afs_family_composite:
    input:
        script = PROJECT_ROOT / "scripts/figures/plot_introner_afs_family_composite.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        helper = PROJECT_ROOT / "scripts/genotyping/plot_introner_2d_afs.py"
    output:
        file0 = FIGURES_DIR / "genotyping" / "introner_afs_family_composite.png"
    params:
        outdir = FIGURES_DIR / "genotyping"
    threads: 1
    log:
        RESULTS / "logs/promoted_analyses/introner_afs_family_composite.log"
    conda: "../envs/analysis_tools.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --genotype-matrix {input.matrix:q} --output {output.file0:q} > {log:q} 2>&1
        """

rule snp_popgen_composite:
    input:
        script = PROJECT_ROOT / "scripts/figures/plot_snp_popgen_composite.py",
        vcf = SNP_DIR / "vcf/mpusilla.snps.4d.notMT.vcf.gz",
        tree = SNP_DIR / "phylogenetics/mpusilla.snps.4d.rooted.treefile",
        model = SNP_DIR / "demography/model_fits.4d.txt",
        bed = SNP_DIR / "degenotate/degeneracy-all-sites.4d.bed.gz",
        colors = PROJECT_ROOT / "master_figure_color_guide.tsv",
        helper = PROJECT_ROOT / "scripts/figure_color_guide.py"
    output:
        file0 = FIGURES_DIR / "snp_popgen" / "snp_popgen_composite.png"
    params:
        outdir = FIGURES_DIR / "snp_popgen"
    threads: 1
    log:
        RESULTS / "logs/promoted_analyses/snp_popgen_composite.log"
    conda: "../envs/analysis_tools.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --vcf {input.vcf:q} --tree {input.tree:q} --model-fits {input.model:q} --target-bed {input.bed:q} --output {output.file0:q} > {log:q} 2>&1
        """

rule figure4_afs_dxy_composite:
    input:
        script = PROJECT_ROOT / "scripts/figures/plot_figure4_afs_dxy_composite.py",
        afs = SNP_DIR / "sfs/afs_by_class.tsv",
        flanks = EVOLUTION_DIR / "diversity_metrics/all_samples_diversity_metrics_200bp.tsv",
        body = EVOLUTION_DIR / "diversity_metrics/shared_introner_body_dxy_200bp.tsv",
        decay = EVOLUTION_DIR / "introner_body_decay/introner_body_decay.per_locus.tsv",
        colors = PROJECT_ROOT / "master_figure_color_guide.tsv",
        helper = PROJECT_ROOT / "scripts/evolution/plot_all_samples_analysis_boxplot.py",
        color_script = PROJECT_ROOT / "scripts/figure_color_guide.py"
    output:
        file0 = FIGURES_DIR / "snp_popgen" / "figure4_afs_dxy_composite.png",
        file1 = FIGURES_DIR / "snp_popgen" / "figure4_population1_pi_supplement.png"
    params:
        outdir = FIGURES_DIR / "snp_popgen"
    threads: 1
    log:
        RESULTS / "logs/promoted_analyses/figure4_afs_dxy_composite.log"
    conda: "../envs/analysis_tools.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --afs {input.afs:q} --flanking-diversity {input.flanks:q} --body-dxy {input.body:q} --body-decay {input.decay:q} --main-output {output.file0:q} --pi-output {output.file1:q} > {log:q} 2>&1
        """

