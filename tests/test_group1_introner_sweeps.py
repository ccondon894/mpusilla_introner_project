from __future__ import annotations

import gzip
import importlib.util
import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TMP_ROOT = Path(os.environ.get("TMPDIR", "/scratch1/chris/tmp"))


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


targets_mod = load_module(
    "build_group1_introner_sweep_targets",
    PROJECT_ROOT / "scripts/popgen/selection_analysis/build_group1_introner_sweep_targets.py",
)
sweeps_mod = load_module(
    "calculate_group1_introner_sweep_windows",
    PROJECT_ROOT / "scripts/popgen/selection_analysis/calculate_group1_introner_sweep_windows.py",
)


@pytest.fixture
def workdir():
    TMP_ROOT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=TMP_ROOT) as path:
        yield Path(path)


def write_matrix(path: Path) -> None:
    group1 = ["CCMP1545", "RCC114", "RCC1614"]
    group2 = ["RCC1749", "RCC3052"]
    rows = []

    def add_ortholog(oid, g1_states, g2_states, status="consistent", contig="CCMP1545#0#scaffold_1", start=1000):
        for sample, presence in zip(group1, g1_states):
            rows.append(
                {
                    "ortholog_id": oid,
                    "sample": sample,
                    "presence": presence,
                    "contig": contig,
                    "start": start,
                    "end": start + 99,
                    "within_group_status": status,
                }
            )
        for sample, presence in zip(group2, g2_states):
            rows.append(
                {
                    "ortholog_id": oid,
                    "sample": sample,
                    "presence": presence,
                    "contig": contig,
                    "start": start,
                    "end": start + 99,
                    "within_group_status": status,
                }
            )

    add_ortholog("keep", [1, 2, 2], [2, 3], start=1000)
    add_ortholog("missing_g1", [1, 2, 3], [2, 2], start=2000)
    add_ortholog("group2_present", [1, 2, 2], [1, 2], start=3000)
    add_ortholog("bad_status", [1, 2, 2], [2, 2], status="low_identity", start=4000)
    add_ortholog("mating_type", [1, 2, 2], [2, 2], contig="CCMP1545#0#scaffold_2", start=1200)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_build_targets_filters_group1_only_and_mating_type(workdir):
    matrix = workdir / "matrix.tsv"
    write_matrix(matrix)
    targets, counts = targets_mod.build_targets(
        genotype_matrix=matrix,
        group1_samples=["CCMP1545", "RCC114", "RCC1614"],
        group2_samples=["RCC1749", "RCC3052"],
        reference_sample="CCMP1545",
        target_mode="group1_only_polymorphic",
        accepted_within_statuses={"consistent", "singleton"},
        mating_contig="CCMP1545#0#scaffold_2",
        mating_start=1000,
        mating_end=2000,
    )
    assert targets["ortholog_id"].tolist() == ["keep"]
    assert counts["group1_complete_target_class_accepted"] == 3
    assert counts["excluded_group2_present"] == 1
    assert counts["excluded_mating_type"] == 1


def test_vcf_parser_tolerates_malformed_format_and_counts_sites(workdir):
    vcf = workdir / "mini.vcf.gz"
    text = "\n".join(
        [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC",
            "ctg\t10\t.\tA\t.\t.\t.\t.\tGT:DP\t0:100,0:100:extra\t0:5\t0:9",
            "ctg\t20\t.\tA\tG\t.\t.\t.\tGT:DP\t1:7\t0:8\t0:9",
            "ctg\t30\t.\tA\tG,C\t.\t.\t.\tGT:DP\t1:7\t0:8\t0:9",
            "ctg\t40\t.\tA\tT\t.\t.\t.\tGT:DP\t.:0\t0:8\t1:9",
            "",
        ]
    )
    with gzip.open(vcf, "wt") as handle:
        handle.write(text)

    sites, stats = sweeps_mod.stream_vcf_sites(vcf, ["A", "B", "C"])
    assert stats["records_seen"] == 4
    assert stats["records_used"] == 2
    assert stats["skipped_multiallelic"] == 1
    assert stats["skipped_missing"] == 1
    assert sites["ctg"]["pos"].tolist() == [10, 20]
    assert sites["ctg"]["alt_count"].tolist() == [0, 1]

    window = sweeps_mod.compute_window_stats(sites["ctg"], "ctg", 1, 25, 3, 1)
    assert window["callable_sites"] == 2
    assert window["segregating_sites"] == 1
    assert window["pass_min_callable"] is True
