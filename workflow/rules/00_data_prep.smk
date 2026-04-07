# ============================================================
# 00_data_prep.smk - Data Preparation Rules
# ============================================================
#
# Prepares input data for downstream analysis:
# - Index reference genomes
# - Create genome dictionaries
# - Validate input files
# - Copy/link assembly files to results
#
# ============================================================

import os
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Paths from config
ASSEMBLIES_SOURCE = Path(config["paths"]["assemblies"])
CCMP1545_REF = Path(config["paths"]["references"]["ccmp1545"])
RCC1749_REF = Path(config["paths"]["references"]["rcc1749"])

# ============================================================
# RULES
# ============================================================

rule copy_assemblies:
    """
    Copy or link sample assemblies (vg_paths.fa) to results directory
    """
    input:
        fa = ASSEMBLIES_SOURCE / "{sample}.vg_paths.fa"
    output:
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    shell:
        """
        cp {input.fa} {output.fa}
        """


rule index_assembly_fasta:
    """
    Index assembly FASTA with samtools faidx
    """
    input:
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        fai = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.fai"
    shell:
        """
        samtools faidx {input.fa}
        """


rule index_assembly_bwa:
    """
    Create BWA index for assembly
    """
    input:
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        bwt = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.bwt",
        pac = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.pac",
        ann = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.ann",
        amb = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.amb",
        sa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.sa"
    threads: 4
    shell:
        """
        bwa index {input.fa}
        """


rule index_assembly_blast:
    """
    Create BLAST database for assembly
    """
    input:
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        nhr = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.nhr",
        nin = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.nin",
        nsq = ASSEMBLIES_DIR / "{sample}.vg_paths.fa.nsq"
    params:
        db_prefix = lambda wildcards, output: str(output.nhr).replace(".nhr", "")
    shell:
        """
        makeblastdb -in {input.fa} -dbtype nucl -out {params.db_prefix}
        """


rule copy_reference_genome:
    """
    Copy reference genome to results
    """
    input:
        ccmp1545 = CCMP1545_REF,
        rcc1749 = RCC1749_REF
    output:
        ccmp1545 = ASSEMBLIES_DIR / "references" / "CCMP1545_v3.fa",
        rcc1749 = ASSEMBLIES_DIR / "references" / "RCC1749_assembly.fa"
    shell:
        """
        mkdir -p $(dirname {output.ccmp1545})
        cp {input.ccmp1545} {output.ccmp1545}
        cp {input.rcc1749} {output.rcc1749}
        """


rule index_reference_bwa:
    """
    Create BWA index for reference genomes
    """
    input:
        fa = ASSEMBLIES_DIR / "references" / "{ref}.fa"
    output:
        bwt = ASSEMBLIES_DIR / "references" / "{ref}.fa.bwt"
    shell:
        """
        bwa index {input.fa}
        """


rule all_data_prep:
    """
    Target: Prepare all data (assemblies, indices)
    """
    input:
        # All sample assemblies copied
        expand(ASSEMBLIES_DIR / "{sample}.vg_paths.fa", sample=ALL_SAMPLES),
        # All FASTA indices
        expand(ASSEMBLIES_DIR / "{sample}.vg_paths.fa.fai", sample=ALL_SAMPLES),
        # All BWA indices
        expand(ASSEMBLIES_DIR / "{sample}.vg_paths.fa.bwt", sample=ALL_SAMPLES),
        # All BLAST databases
        expand(ASSEMBLIES_DIR / "{sample}.vg_paths.fa.nhr", sample=ALL_SAMPLES),
        # Reference genomes
        ASSEMBLIES_DIR / "references" / "CCMP1545_v3.fa",
        ASSEMBLIES_DIR / "references" / "RCC1749_assembly.fa"


rule validate_inputs:
    """
    Validate that all required input files exist
    """
    output:
        touch(RESULTS / ".inputs_validated")
    run:
        import sys

        errors = []

        # Check assemblies
        for sample in ALL_SAMPLES:
            fa = ASSEMBLIES_SOURCE / f"{sample}.vg_paths.fa"
            if not fa.exists():
                errors.append(f"Missing assembly: {fa}")

        # Check references
        if not Path(CCMP1545_REF).exists():
            errors.append(f"Missing CCMP1545 reference: {CCMP1545_REF}")
        if not Path(RCC1749_REF).exists():
            errors.append(f"Missing RCC1749 reference: {RCC1749_REF}")

        # Check introner truth set
        truth_set = Path(PROJECT_ROOT) / "resources" / "introner_truth_set.fa"
        if not truth_set.exists():
            errors.append(f"Missing introner truth set: {truth_set}")

        if errors:
            print("Validation failed:")
            for e in errors:
                print(f"  - {e}")
            sys.exit(1)
        else:
            print("All inputs validated successfully!")
