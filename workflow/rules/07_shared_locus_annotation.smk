# Manuscript-relevant analysis promoted from analysis/.

rule unresolved_shared_locus_annotation:
    input:
        script = PROJECT_ROOT / "scripts/genotyping/qa/summarize_fingerprint_unresolved_shared_loci.py",
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        fingerprints = expand(GENOTYPING_DIR / "insertion_fingerprints/{sample}.fingerprints.tsv", sample=ALL_SAMPLES),
        gtfs = expand(ANNOTATIONS_DIR / "{sample}.gtf", sample=ALL_SAMPLES)
    output:
        report = directory(GENOTYPING_DIR / "qa/unresolved_shared_loci")
    params:
        outdir = GENOTYPING_DIR / "qa/unresolved_shared_loci"
    threads: 1
    log:
        RESULTS / "logs/promoted_analyses/unresolved_shared_locus_annotation.log"
    conda: "../envs/analysis_tools.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} $(dirname {log:q})
        export TMPDIR=/scratch1/chris/tmp MPLCONFIGDIR=/scratch1/chris/tmp/matplotlib
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
        python {input.script:q} --matrix {input.matrix:q} --fingerprint-dir {GENOTYPING_DIR:q}/insertion_fingerprints --annotation-dir {ANNOTATIONS_DIR:q} --output-prefix {params.outdir:q}/fingerprint_unresolved_shared_loci > {log:q} 2>&1
        """

