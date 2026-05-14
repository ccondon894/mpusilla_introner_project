# ============================================================
# 30_expression.smk - RNA-seq Expression Quantification
# ============================================================
#
# Performs RNA-seq alignment and expression quantification:
# 1. HISAT2 alignment of RNA-seq reads
# 2. featureCounts gene-level quantification
# 3. StringTie de novo transcript assembly & GTF augmentation
# 4. SQANTI3 long-read isoform classification (R2C2 data)
#
# Adapted from:
# - /scratch1/chris/introner-expression-analysis/rules/hisat_2.smk
# - /scratch1/chris/introner-expression-analysis/rules/create_count_matrix.smk
#
# ============================================================

import os
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Output directories
EXPRESSION_DIR = RESULTS / "expression"
ALIGNMENT_DIR = EXPRESSION_DIR / "alignments"
COUNTS_DIR = EXPRESSION_DIR / "counts"
SQANTI_DIR = EXPRESSION_DIR / "sqanti3"
STRINGTIE_DIR = EXPRESSION_DIR / "stringtie"
EXPRESSION_LOG_DIR = EXPRESSION_DIR / "logs"

# Input paths
RNA_SEQ_DIR = Path(config["paths"]["original"]["rna_seq"])
R2C2_BAM_DIR = Path(config["paths"]["original"]["r2c2_bams"])

# Samples with RNA-seq data
RNA_SEQ_SAMPLES = config["samples"]["with_rna_seq"]

# RNA-seq replicates per sample
REPLICATES = {
    "CCMP1545": ["834_A1", "834_A2", "834_A3", "834_B1"],
    "RCC1614": ["1614_A1", "1614_A2", "1614_A3", "1614_A4"],
    "RCC1749": ["1749_A1", "1749_A2", "1749_A4", "1749_B1"]
}


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_all_replicates():
    """Get flat list of all replicate names"""
    replicates = []
    for sample, reps in REPLICATES.items():
        replicates.extend(reps)
    return replicates


def get_sample_from_replicate(replicate):
    """Map replicate name to sample"""
    for sample, reps in REPLICATES.items():
        if replicate in reps:
            return sample
    return None


# ============================================================
# HISAT2 INDEXING (per-sample)
# ============================================================

rule hisat2_index:
    """
    Build HISAT2 index for each sample's VG paths assembly.

    Each strain is aligned to its own assembly so contig names
    match the lifted-over GTF annotations.
    """
    input:
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        index = directory(ALIGNMENT_DIR / "hisat2_index_{sample}"),
        done = ALIGNMENT_DIR / "hisat2_index_{sample}" / "index.done"
    params:
        prefix = lambda wildcards: str(ALIGNMENT_DIR / f"hisat2_index_{wildcards.sample}" / wildcards.sample)
    threads: 8
    log:
        EXPRESSION_LOG_DIR / "hisat2_index_{sample}.log"
    conda: "../envs/expression.yaml"
    shell:
        """
        mkdir -p {output.index}
        mkdir -p {EXPRESSION_LOG_DIR}
        hisat2-build -p {threads} {input.fa} {params.prefix} 2> {log}
        touch {output.done}
        """


# ============================================================
# RNA-SEQ ALIGNMENT (per-sample index)
# ============================================================

rule hisat2_align:
    """
    Align RNA-seq reads with HISAT2.

    Each replicate is aligned to its parent strain's VG paths assembly.
    """
    input:
        index_done = lambda wildcards: ALIGNMENT_DIR / f"hisat2_index_{get_sample_from_replicate(wildcards.replicate)}" / "index.done",
        r1 = lambda wildcards: RNA_SEQ_DIR / get_sample_from_replicate(wildcards.replicate) / f"{wildcards.replicate}_R1.fastq.gz",
        r2 = lambda wildcards: RNA_SEQ_DIR / get_sample_from_replicate(wildcards.replicate) / f"{wildcards.replicate}_R2.fastq.gz"
    output:
        bam = ALIGNMENT_DIR / "bams" / "{replicate}.sorted.bam",
        bai = ALIGNMENT_DIR / "bams" / "{replicate}.sorted.bam.bai"
    params:
        index_prefix = lambda wildcards: str(ALIGNMENT_DIR / f"hisat2_index_{get_sample_from_replicate(wildcards.replicate)}" / get_sample_from_replicate(wildcards.replicate))
    threads: 8
    log:
        EXPRESSION_LOG_DIR / "hisat2_align" / "{replicate}.log"
    conda: "../envs/expression.yaml"
    shell:
        """
        mkdir -p {ALIGNMENT_DIR}/bams
        mkdir -p {EXPRESSION_LOG_DIR}/hisat2_align

        hisat2 -p {threads} \
            -x {params.index_prefix} \
            -1 {input.r1} \
            -2 {input.r2} \
            2> {log} | \
        samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam}
        """


# ============================================================
# FEATURECOUNTS QUANTIFICATION
# ============================================================

rule featurecounts_per_sample:
    """
    Run featureCounts on individual BAM files.

    Counts reads mapping to exons, summarized at gene level.
    Uses each strain's own GTF to match VG paths contig names.
    """
    input:
        bam = ALIGNMENT_DIR / "bams" / "{replicate}.sorted.bam",
        gtf = lambda wildcards: ANNOTATIONS_DIR / f"{get_sample_from_replicate(wildcards.replicate)}.gtf"
    output:
        counts = COUNTS_DIR / "{replicate}.counts.txt",
        summary = COUNTS_DIR / "{replicate}.counts.txt.summary"
    threads: 4
    log:
        EXPRESSION_LOG_DIR / "featurecounts" / "{replicate}.log"
    conda: "../envs/expression.yaml"
    shell:
        """
        mkdir -p {COUNTS_DIR}
        mkdir -p {EXPRESSION_LOG_DIR}/featurecounts

        featureCounts \
            -a {input.gtf} \
            -o {output.counts} \
            -t exon \
            -g gene_id \
            -T {threads} \
            -p \
            {input.bam} \
            2> {log}
        """


rule merge_count_matrix:
    """
    Merge individual count files into a single count matrix.
    """
    input:
        counts = expand(COUNTS_DIR / "{replicate}.counts.txt", replicate=get_all_replicates())
    output:
        matrix = COUNTS_DIR / "merged_counts_matrix.csv"
    run:
        import pandas as pd

        # Read first file to get gene IDs
        # Note: use skiprows=1 instead of comment='#' because contig
        # names contain '#' (e.g. CCMP1545#0#scaffold_1) which causes
        # comment='#' to corrupt parsing
        first_file = input.counts[0]
        df = pd.read_csv(first_file, sep='\t', skiprows=1, index_col=0)

        # Keep only Geneid and count column
        result = df[['Length']].copy()

        # Add counts from each file
        for count_file in input.counts:
            replicate = os.path.basename(count_file).replace('.counts.txt', '')
            temp_df = pd.read_csv(count_file, sep='\t', skiprows=1, index_col=0)
            # Count column is the last one (BAM filename)
            count_col = temp_df.columns[-1]
            result[replicate] = temp_df[count_col]

        result.to_csv(output.matrix)


# ============================================================
# STRINGTIE DE NOVO ANNOTATION
# ============================================================
# Augments lifted-over GTF annotations with de novo transcript
# models assembled from RNA-seq data. This recovers genes that
# were missed during miniprot-based liftover due to sequence
# divergence or lineage-specific genes.
#
# Pipeline:
# 1. Per-replicate guided assembly (stringtie_assemble)
# 2. Merge assemblies per sample (stringtie_merge)
# 3. Filter & augment the existing GTF (augment_gtf)
# ============================================================

rule stringtie_assemble:
    """
    Reference-guided transcript assembly per replicate.

    Uses the existing lifted-over GTF as a guide (-G) so known
    gene structures are preserved, while also assembling novel
    transcripts from RNA-seq coverage in unannotated regions.
    """
    input:
        bam = ALIGNMENT_DIR / "bams" / "{replicate}.sorted.bam",
        gtf = lambda wildcards: ANNOTATIONS_DIR / f"{get_sample_from_replicate(wildcards.replicate)}.gtf"
    output:
        gtf = STRINGTIE_DIR / "per_replicate" / "{replicate}.stringtie.gtf"
    threads: 8
    log:
        EXPRESSION_LOG_DIR / "stringtie" / "{replicate}.log"
    conda: "../envs/expression.yaml"
    shell:
        """
        mkdir -p $(dirname {output.gtf})
        mkdir -p $(dirname {log})

        stringtie {input.bam} \
            -G {input.gtf} \
            -o {output.gtf} \
            -p {threads} \
            -f 0.05 \
            -m 200 \
            -c 2.5 \
            2> {log}
        """


rule stringtie_merge:
    """
    Merge per-replicate StringTie assemblies for each sample.

    Combines transcript models across replicates to produce a
    unified annotation. The existing GTF is included as a guide
    so all original annotations are retained in the merged output.
    """
    input:
        rep_gtfs = lambda wildcards: expand(
            STRINGTIE_DIR / "per_replicate" / "{replicate}.stringtie.gtf",
            replicate=REPLICATES[wildcards.sample]
        ),
        guide_gtf = ANNOTATIONS_DIR / "{sample}.gtf"
    output:
        gtf = STRINGTIE_DIR / "merged" / "{sample}.stringtie_merged.gtf"
    threads: 8
    log:
        EXPRESSION_LOG_DIR / "stringtie" / "{sample}_merge.log"
    conda: "../envs/expression.yaml"
    shell:
        """
        mkdir -p $(dirname {output.gtf})

        stringtie --merge \
            -G {input.guide_gtf} \
            -o {output.gtf} \
            -T 1.0 \
            -f 0.05 \
            -i \
            {input.rep_gtfs} \
            2> {log}
        """


rule augment_gtf:
    """
    Filter novel StringTie transcripts and merge into existing GTF.

    Identifies transcripts from StringTie that do NOT overlap any
    existing annotation (novel genes), filters them by minimum
    expression and structure, then appends them to the original GTF
    to produce an augmented annotation.
    """
    input:
        merged_gtf = STRINGTIE_DIR / "merged" / "{sample}.stringtie_merged.gtf",
        original_gtf = ANNOTATIONS_DIR / "{sample}.gtf"
    output:
        augmented_gtf = ANNOTATIONS_DIR / "{sample}.augmented.gtf",
        novel_gtf = STRINGTIE_DIR / "novel" / "{sample}.novel_transcripts.gtf",
        stats = STRINGTIE_DIR / "novel" / "{sample}.augmentation_stats.txt"
    params:
        min_length = 200,
        min_exons = 1
    log:
        EXPRESSION_LOG_DIR / "stringtie" / "{sample}_augment.log"
    run:
        import re
        from collections import defaultdict

        def parse_attribute(line, key):
            """Extract an attribute value from a GTF line."""
            match = re.search(rf'{key} "([^"]*)"', line)
            return match.group(1) if match else None

        def parse_gtf_intervals(gtf_path):
            """Parse GTF into per-chromosome interval list."""
            intervals = defaultdict(list)
            with open(gtf_path) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) < 9:
                        continue
                    if fields[2] == "transcript":
                        chrom = fields[0]
                        start, end = int(fields[3]), int(fields[4])
                        intervals[chrom].append((start, end))
            # Sort intervals for binary search
            for chrom in intervals:
                intervals[chrom].sort()
            return intervals

        def overlaps_any(chrom, start, end, intervals):
            """Check if a region overlaps any interval on the chromosome."""
            if chrom not in intervals:
                return False
            import bisect
            ivs = intervals[chrom]
            # Find intervals that could overlap
            idx = bisect.bisect_right(ivs, (end,)) - 1
            # Check nearby intervals
            for i in range(max(0, idx - 5), min(len(ivs), idx + 5)):
                iv_start, iv_end = ivs[i]
                if iv_start <= end and iv_end >= start:
                    return True
            return False

        # Parse existing annotations
        original_intervals = parse_gtf_intervals(str(input.original_gtf))

        # Parse merged StringTie output, group lines by transcript
        transcripts = defaultdict(list)
        transcript_info = {}
        with open(str(input.merged_gtf)) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue
                tid = parse_attribute(line, "transcript_id")
                if tid:
                    transcripts[tid].append(line)
                    if fields[2] == "transcript":
                        chrom = fields[0]
                        start, end = int(fields[3]), int(fields[4])
                        transcript_info[tid] = {
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "strand": fields[6],
                            "length": end - start + 1
                        }
                    elif fields[2] == "exon":
                        transcript_info.setdefault(tid, {})
                        transcript_info[tid].setdefault("n_exons", 0)
                        transcript_info[tid]["n_exons"] += 1

        # Identify novel transcripts (no overlap with original GTF)
        novel_tids = []
        for tid, info in transcript_info.items():
            if "chrom" not in info:
                continue
            if info["length"] < int(params.min_length):
                continue
            if info.get("n_exons", 1) < int(params.min_exons):
                continue
            if not overlaps_any(info["chrom"], info["start"], info["end"],
                                original_intervals):
                novel_tids.append(tid)

        # Write novel transcripts GTF
        os.makedirs(os.path.dirname(str(output.novel_gtf)), exist_ok=True)
        novel_gene_ids = set()
        with open(str(output.novel_gtf), "w") as f:
            for tid in sorted(novel_tids):
                for line in transcripts[tid]:
                    # Prefix novel gene_ids with "NOVEL_" for clarity
                    gid = parse_attribute(line, "gene_id")
                    if gid:
                        novel_gene_ids.add(gid)
                    f.write(line)

        # Write augmented GTF: original + novel
        with open(str(output.augmented_gtf), "w") as out:
            # Copy original GTF
            with open(str(input.original_gtf)) as orig:
                for line in orig:
                    out.write(line)
            # Append novel transcripts
            out.write(f"\n# Novel transcripts from StringTie de novo assembly\n")
            with open(str(output.novel_gtf)) as novel:
                for line in novel:
                    out.write(line)

        # Write stats
        original_genes = set()
        with open(str(input.original_gtf)) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                gid = parse_attribute(line, "gene_id")
                if gid:
                    original_genes.add(gid)

        with open(str(output.stats), "w") as f:
            f.write(f"Original genes: {len(original_genes)}\n")
            f.write(f"Novel transcripts added: {len(novel_tids)}\n")
            f.write(f"Novel gene loci added: {len(novel_gene_ids)}\n")
            f.write(f"Total genes in augmented GTF: "
                    f"{len(original_genes) + len(novel_gene_ids)}\n")
            f.write(f"\nFilters applied:\n")
            f.write(f"  Min transcript length: {params.min_length} bp\n")
            f.write(f"  Min exons: {params.min_exons}\n")


rule find_novel_orthologs:
    """
    BLASTX novel StringTie transcripts against CCMP1545 proteins
    to identify orthologs, then rename gene_ids in the augmented GTF
    to match the CCMP1545 gene names.

    Produces:
    - BLASTX results table
    - Augmented GTF with CCMP1545 gene names for orthologs
    - Ortholog mapping summary
    """
    input:
        novel_gtf = STRINGTIE_DIR / "novel" / "{sample}.novel_transcripts.gtf",
        augmented_gtf = ANNOTATIONS_DIR / "{sample}.augmented.gtf",
        assembly = ASSEMBLIES_DIR / "{sample}.vg_paths.fa",
        proteins = RESULTS / "proteins" / "reference_proteins.fa"
    output:
        blastx = STRINGTIE_DIR / "orthologs" / "{sample}.blastx.tsv"
    params:
        evalue = "1e-5",
        blastdb = lambda wildcards: str(STRINGTIE_DIR / "orthologs" / f"{wildcards.sample}_ccmp1545_proteins")
    threads: 8
    log:
        EXPRESSION_LOG_DIR / "stringtie" / "{sample}_orthologs.log"
    conda: "../envs/expression.yaml"
    shell:
        """
        mkdir -p $(dirname {output.blastx})
        mkdir -p $(dirname {log})

        # Extract novel transcript sequences
        gffread {input.novel_gtf} -g {input.assembly} \
            -w {output.blastx}.tmp.fa 2>> {log}

        # Build BLAST database
        makeblastdb -in {input.proteins} -dbtype prot \
            -out {params.blastdb} 2>> {log}

        # Run BLASTX
        blastx -query {output.blastx}.tmp.fa \
            -db {params.blastdb} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
            -evalue {params.evalue} \
            -max_target_seqs 1 \
            -num_threads {threads} \
            -out {output.blastx} 2>> {log}

        rm -f {output.blastx}.tmp.fa

        # Clean up BLAST database files
        rm -f {params.blastdb}.pdb {params.blastdb}.phr {params.blastdb}.pin \
              {params.blastdb}.pjs {params.blastdb}.pot {params.blastdb}.psq \
              {params.blastdb}.ptf {params.blastdb}.pto
        """


rule rename_novel_genes:
    """
    Rename novel gene_ids in the augmented GTF to match CCMP1545
    orthologs identified by BLASTX.

    For each novel transcript with a BLASTX hit, the gene_id is
    updated to the CCMP1545 gene name (with _novel suffix to
    distinguish from the miniprot-derived annotations). Transcripts
    with no ortholog retain their StringTie MSTRG.* identifiers.
    """
    input:
        blastx = STRINGTIE_DIR / "orthologs" / "{sample}.blastx.tsv",
        augmented_gtf = ANNOTATIONS_DIR / "{sample}.augmented.gtf",
        proteins = RESULTS / "proteins" / "reference_proteins.fa"
    output:
        renamed_gtf = ANNOTATIONS_DIR / "{sample}.augmented.renamed.gtf",
        mapping = STRINGTIE_DIR / "orthologs" / "{sample}.ortholog_mapping.tsv",
        stats = STRINGTIE_DIR / "orthologs" / "{sample}.ortholog_stats.txt"
    run:
        import re
        from collections import defaultdict

        # Build transcript_id -> gene_id mapping from protein FASTA headers
        # Headers look like: >58671.3.0.228 gene=MicpuC2.EuGene.0000060180.3.0.228 ...
        prot_to_gene = {}
        with open(str(input.proteins)) as f:
            for line in f:
                if line.startswith(">"):
                    prot_id = line.split()[0][1:]  # remove >
                    gene_match = re.search(r'gene=(\S+)', line)
                    if gene_match:
                        prot_to_gene[prot_id] = gene_match.group(1)

        # Parse BLASTX results: pick best hit per query transcript
        # (highest bitscore)
        best_hits = {}
        with open(str(input.blastx)) as f:
            for line in f:
                fields = line.strip().split("\t")
                qid = fields[0]       # StringTie transcript_id
                sid = fields[1]       # CCMP1545 protein/transcript_id
                pident = float(fields[2])
                bitscore = float(fields[11])

                if qid not in best_hits or bitscore > best_hits[qid][2]:
                    best_hits[qid] = (sid, pident, bitscore)

        # Map StringTie transcript_id -> CCMP1545 gene_id
        # StringTie GTF transcript_ids include the gene_id (MSTRG.X) as part
        # of the transcript_id (MSTRG.X.Y), so we need to map at both levels
        tid_to_ccmp_gene = {}
        for qid, (sid, pident, bitscore) in best_hits.items():
            ccmp_gene = prot_to_gene.get(sid, sid)
            tid_to_ccmp_gene[qid] = (ccmp_gene, pident, bitscore)

        # Build MSTRG gene_id -> CCMP1545 gene_id mapping
        # A MSTRG gene may have multiple transcripts; use best hit
        mstrg_to_ccmp = {}
        for tid, (ccmp_gene, pident, bitscore) in tid_to_ccmp_gene.items():
            # Extract MSTRG gene_id from transcript_id (MSTRG.123.1 -> MSTRG.123)
            parts = tid.rsplit(".", 1)
            mstrg_gene = parts[0] if len(parts) > 1 else tid

            if mstrg_gene not in mstrg_to_ccmp or bitscore > mstrg_to_ccmp[mstrg_gene][2]:
                mstrg_to_ccmp[mstrg_gene] = (ccmp_gene, pident, bitscore)

        # Write mapping file
        with open(str(output.mapping), "w") as f:
            f.write("stringtie_gene_id\tccmp1545_gene_id\tpident\tbitscore\n")
            for mstrg, (ccmp, pident, bitscore) in sorted(mstrg_to_ccmp.items()):
                f.write(f"{mstrg}\t{ccmp}\t{pident:.1f}\t{bitscore:.1f}\n")

        # Rename gene_ids in augmented GTF
        # Track which CCMP1545 gene names are already in the original part
        # of the GTF to handle potential duplicates
        existing_genes = set()
        renamed_count = 0
        no_ortholog_count = 0
        novel_section = False

        with open(str(output.renamed_gtf), "w") as out:
            with open(str(input.augmented_gtf)) as f:
                for line in f:
                    if "Novel transcripts from StringTie" in line:
                        novel_section = True
                        out.write(line)
                        continue

                    if line.startswith("#") or line.strip() == "":
                        out.write(line)
                        continue

                    if not novel_section:
                        # Original annotation — track gene_ids, pass through
                        gid_match = re.search(r'gene_id "([^"]*)"', line)
                        if gid_match:
                            existing_genes.add(gid_match.group(1))
                        out.write(line)
                    else:
                        # Novel transcript — rename if ortholog found
                        gid_match = re.search(r'gene_id "([^"]*)"', line)
                        if gid_match:
                            mstrg_gene = gid_match.group(1)
                            if mstrg_gene in mstrg_to_ccmp:
                                ccmp_gene = mstrg_to_ccmp[mstrg_gene][0]
                                # Add _novel suffix to distinguish from
                                # miniprot-derived annotation
                                new_gene_id = ccmp_gene + "_novel"
                                line = line.replace(
                                    f'gene_id "{mstrg_gene}"',
                                    f'gene_id "{new_gene_id}"'
                                )
                                renamed_count += 1
                            else:
                                no_ortholog_count += 1
                        out.write(line)

        # Write stats
        with open(str(output.stats), "w") as f:
            f.write(f"BLASTX hits (unique transcripts): {len(best_hits)}\n")
            f.write(f"Novel gene loci with ortholog: {len(mstrg_to_ccmp)}\n")
            f.write(f"Novel gene loci without ortholog: "
                    f"{no_ortholog_count}\n")
            f.write(f"GTF lines renamed: {renamed_count}\n")
            f.write(f"\nIdentity distribution of orthologs:\n")
            brackets = {"<50%": 0, "50-70%": 0, "70-90%": 0, "90-100%": 0}
            for _, (_, pident, _) in mstrg_to_ccmp.items():
                if pident < 50:
                    brackets["<50%"] += 1
                elif pident < 70:
                    brackets["50-70%"] += 1
                elif pident < 90:
                    brackets["70-90%"] += 1
                else:
                    brackets["90-100%"] += 1
            for bracket, count in brackets.items():
                f.write(f"  {bracket}: {count}\n")


# ============================================================
# SQANTI3 LONG-READ ANALYSIS
# ============================================================

rule parse_sqanti3_output:
    """
    Parse SQANTI3 classification files for isoform analysis.

    Extracts isoform classifications, structural categories, and
    abundance information from R2C2 long-read data.
    """
    input:
        # SQANTI3 classification files (pre-computed, copied to data/)
        ccmp1545 = PROJECT_ROOT / "data" / "sqanti3_output_834" / "834_isoforms_classification.filtered.txt",
        rcc1614 = PROJECT_ROOT / "data" / "sqanti3_output_1614" / "1614_isoforms_classification.filtered.txt",
        rcc1749 = PROJECT_ROOT / "data" / "sqanti3_output_1749" / "1749_isoforms_classification.filtered.txt"
    output:
        parsed = SQANTI_DIR / "parsed_sqanti3_data.tsv",
        summary = SQANTI_DIR / "sqanti3_summary.txt"
    log:
        EXPRESSION_LOG_DIR / "parse_sqanti3.log"
    conda: "../envs/expression.yaml"
    shell:
        """
        mkdir -p {SQANTI_DIR}

        python {PROJECT_ROOT}/scripts/expression/sqanti3/parse_sqanti3_data.py \
            --ccmp1545 {input.ccmp1545} \
            --rcc1614 {input.rcc1614} \
            --rcc1749 {input.rcc1749} \
            --output {output.parsed} \
            --summary {output.summary} \
            2> {log}
        """


# ============================================================
# TARGET RULES
# ============================================================

rule expression_complete:
    """
    Target: Complete expression quantification pipeline.
    """
    input:
        COUNTS_DIR / "merged_counts_matrix.csv",
        SQANTI_DIR / "parsed_sqanti3_data.tsv"


rule alignment_only:
    """
    Target: Run HISAT2 alignments only.
    """
    input:
        expand(ALIGNMENT_DIR / "bams" / "{replicate}.sorted.bam", replicate=get_all_replicates())


rule counts_only:
    """
    Target: Generate count matrix only.
    """
    input:
        COUNTS_DIR / "merged_counts_matrix.csv"


rule stringtie_denovo:
    """
    Target: Run StringTie de novo annotation for all RNA-seq samples.

    Produces augmented + renamed GTF files at
    results/annotations/{sample}.augmented.renamed.gtf
    """
    input:
        expand(ANNOTATIONS_DIR / "{sample}.augmented.renamed.gtf", sample=RNA_SEQ_SAMPLES)
