#!/usr/bin/env python3
"""
Consolidated ortholog group resolver (Milestone 3).

Iterates classify -> split -> reclassify until stable, then applies sequence
identity refinement. Replaces the multi-rule chain:
  classify_sharing_status -> split_overmerged_orthologs -> reclassify_after_split
  -> split_cross_group_mispairs_pass2 -> reclassify_after_cross_split
  -> compare_introner_sequences
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import classify_sharing_status as css
import compare_introner_sequences as seqcmp
import split_overmerged_orthologs as splitter

GROUP2_SAMPLES = css.GROUP2_SAMPLES
CROSS_GROUP_SPLIT_REASONS = splitter.CROSS_GROUP_SPLIT_REASONS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resolve ortholog groups in one iterative pass."
    )
    parser.add_argument("--matrix", required=True, help="Oriented genotype matrix TSV")
    parser.add_argument("--fingerprints", required=True, nargs="+")
    parser.add_argument("--beds", required=True, nargs="+")
    parser.add_argument("--gtfs", required=True, nargs="+")
    parser.add_argument("--tandems", nargs="*", default=[])
    parser.add_argument("--genome-dir", required=True)
    parser.add_argument("--locus-group-members", default=None,
                        help="Optional locus_group_members.tsv sidecar from network build")
    parser.add_argument("--output", required=True,
                        help="Resolved pre-annotation genotype matrix TSV")
    parser.add_argument("--insertion-group-map", required=True)
    parser.add_argument("--independent-insertion-events", required=True)
    parser.add_argument("--resolution-summary", required=True)
    parser.add_argument("--codon-tolerance", type=int, default=css.CODON_TOLERANCE)
    parser.add_argument("--flanking-length", type=int, default=100)
    parser.add_argument("--within-group-identity-threshold", type=float, default=0.80)
    parser.add_argument("--cross-group-identity-threshold", type=float, default=0.60)
    parser.add_argument("--max-iterations", type=int, default=10)
    return parser.parse_args()


def read_matrix(path: str) -> tuple[list[dict], list[str]]:
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    return rows, fieldnames


def write_matrix(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def group_rows_by_ortholog(rows: list[dict]) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = defaultdict(list)
    for idx, row in enumerate(rows):
        groups[row["ortholog_id"]].append(idx)
    return groups


def classify_matrix(
    rows: list[dict],
    fieldnames: list[str],
    fps: dict,
    codon_tolerance: int,
) -> tuple[list[dict], list[str], list[dict]]:
    ortholog_groups = group_rows_by_ortholog(rows)
    group_summaries: list[dict] = []

    for oid, row_indices in ortholog_groups.items():
        members = []
        for idx in row_indices:
            row = rows[idx]
            if row["presence"] != "1":
                continue
            sample = row["sample"]
            group = "G2" if sample in GROUP2_SAMPLES else "G1"
            key = (sample, row["contig"], int(row["start"]), int(row["end"]))
            fp = fps.get(key, {})
            members.append(
                {
                    "sample": sample,
                    "group": group,
                    "family": row.get("family", ""),
                    "location_type": fp.get("location_type", ""),
                    "codon_number": fp.get("codon_number", ""),
                    "codon_offset": fp.get("codon_offset", ""),
                    "exon_number": fp.get("exon_number", ""),
                    "intron_number": fp.get("intron_number", ""),
                    "flanking_aa_context": fp.get("flanking_aa_context", ""),
                    "confidence": fp.get("confidence", "no_fingerprint"),
                    "gene_id": fp.get("gene_id", ""),
                }
            )

        within_status, cross_status, cross_reason = css.classify_ortholog_group(
            members, codon_tolerance
        )

        if members:
            typed = [m for m in members if m["confidence"] == "high"]

            def format_locus(member: dict) -> str:
                ctx = member.get("flanking_aa_context", "")
                if ctx:
                    return f"{member['sample']}:{ctx}"
                if member["location_type"] == "cds":
                    return (
                        f"{member['sample']}:cds.{member['codon_number']}."
                        f"{member['codon_offset']}"
                    )
                if member["location_type"] == "intron":
                    return f"{member['sample']}:intron.{member['intron_number']}"
                return f"{member['sample']}:?"

            group_summaries.append(
                {
                    "ortholog_id": oid,
                    "within_group_status": within_status,
                    "cross_group_status": cross_status,
                    "cross_group_reason": cross_reason,
                    "n_present": len(members),
                    "n_fingerprinted": len(typed),
                    "families": ";".join(sorted({m["family"] for m in members})),
                    "loci": ";".join(format_locus(m) for m in typed),
                }
            )

        for idx in row_indices:
            rows[idx]["within_group_status"] = within_status
            rows[idx]["cross_group_status"] = cross_status
            rows[idx]["cross_group_reason"] = cross_reason

    out_fieldnames = [
        f
        for f in fieldnames
        if f
        not in (
            "sharing_status",
            "within_group_status",
            "cross_group_status",
            "cross_group_reason",
        )
    ]
    out_fieldnames += [
        "within_group_status",
        "cross_group_status",
        "cross_group_reason",
    ]
    return rows, out_fieldnames, group_summaries


def split_matrix(
    rows: list[dict],
    fieldnames: list[str],
    fps_by_coord: dict,
    fps_by_sample_gene: dict,
    beds_by_coord: dict,
    tandems: dict,
    cds_by_sample_gene: dict,
    codon_tolerance: int,
    cross_group_only: bool,
    iteration: int,
) -> tuple[list[dict], list[str], list[dict], int]:
    ortholog_groups = group_rows_by_ortholog(rows)
    output_rows: list[dict] = []
    mapping_records: list[dict] = []
    n_split = 0

    for oid in sorted(ortholog_groups):
        row_indices = ortholog_groups[oid]
        if cross_group_only:
            cross_reason = rows[row_indices[0]].get("cross_group_reason", "")
            if cross_reason in CROSS_GROUP_SPLIT_REASONS:
                result = splitter.split_cross_group_mispair(
                    oid,
                    row_indices,
                    rows,
                    fps_by_coord,
                    cds_by_sample_gene,
                    codon_tolerance,
                )
            else:
                result = splitter.split_group_by_existing_family(
                    oid, row_indices, rows
                )
        else:
            result = splitter.split_group_with_recovery(
                oid,
                row_indices,
                rows,
                fps_by_coord,
                fps_by_sample_gene,
                beds_by_coord,
                tandems,
                cds_by_sample_gene,
                codon_tolerance,
            )

        if result is None:
            for idx in row_indices:
                output_rows.append(rows[idx])
            n_present = sum(1 for idx in row_indices if rows[idx]["presence"] == "1")
            mapping_records.append(
                {
                    "iteration": iteration,
                    "resolution_pass": "cross_group" if cross_group_only else "within_group",
                    "original_id": oid,
                    "new_id": oid,
                    "n_clusters": 1,
                    "n_members": n_present,
                    "n_recovered": 0,
                    "has_tandem": "",
                    "note": "unchanged",
                }
            )
            continue

        n_split += 1
        n_clusters = len(result)
        original_pres1 = sum(
            1 for idx in row_indices if rows[idx]["presence"] == "1"
        )
        new_pres1 = sum(
            1 for _, split_rows, _, _ in result for r in split_rows if r["presence"] == "1"
        )
        recovered_in_group = new_pres1 - original_pres1

        if any("_cgsplit_" in item[0] for item in result):
            note = "cross_group_split"
        elif any("_famsplit" in item[0] for item in result):
            note = "family_split"
        else:
            note = "split"

        for new_id, split_rows, _, has_tandem in result:
            output_rows.extend(split_rows)
            mapping_records.append(
                {
                    "iteration": iteration,
                    "resolution_pass": "cross_group" if cross_group_only else "within_group",
                    "original_id": oid,
                    "new_id": new_id,
                    "n_clusters": n_clusters,
                    "n_members": sum(1 for r in split_rows if r["presence"] == "1"),
                    "n_recovered": recovered_in_group if new_id == result[0][0] else "",
                    "has_tandem": "1" if has_tandem else "0",
                    "note": note,
                }
            )

    out_fieldnames = list(fieldnames)
    if "original_ortholog_id" not in out_fieldnames:
        out_fieldnames.append("original_ortholog_id")
    for row in output_rows:
        row.setdefault("original_ortholog_id", row.get("ortholog_id", ""))

    return output_rows, out_fieldnames, mapping_records, n_split


def compare_matrix(
    rows: list[dict],
    fieldnames: list[str],
    genome_dir: str,
    flanking_length: int,
    within_threshold: float,
    cross_threshold: float,
) -> tuple[list[dict], list[str], list[dict]]:
    ortholog_groups = group_rows_by_ortholog(rows)
    all_samples = sorted({row["sample"] for row in rows})
    genomes = seqcmp.load_genomes(genome_dir, all_samples)
    aligner = seqcmp.make_aligner()
    group_summaries: list[dict] = []

    try:
        for oid, row_indices in ortholog_groups.items():
            first_row = rows[row_indices[0]]
            within_status = first_row.get("within_group_status", "")
            cross_status = first_row.get("cross_group_status", "")

            members = []
            for idx in row_indices:
                row = rows[idx]
                if row["presence"] != "1":
                    continue
                sample = row["sample"]
                members.append(
                    {
                        "sample": sample,
                        "contig": row["contig"],
                        "start": int(row["start"]),
                        "end": int(row["end"]),
                        "group": "G2" if sample in GROUP2_SAMPLES else "G1",
                        "family": row.get("family", ""),
                    }
                )

            within_ident = None
            cross_ident = None
            n_extracted = 0
            if len(members) >= 2:
                within_ident, cross_ident, n_extracted = seqcmp.compute_group_identities(
                    members, genomes, aligner, flanking_length
                )

            families_consistent = seqcmp.check_family_consistency(members)
            refined_within = seqcmp.refine_within_group(
                within_status,
                within_ident,
                within_threshold,
                families_consistent,
            )
            refined_cross = seqcmp.refine_cross_group(
                cross_status, cross_ident, cross_threshold
            )

            cross_str = f"{cross_ident:.4f}" if cross_ident is not None else ""
            within_str = f"{within_ident:.4f}" if within_ident is not None else ""

            for idx in row_indices:
                rows[idx]["within_group_status"] = refined_within
                rows[idx]["cross_group_status"] = refined_cross
                rows[idx]["within_group_identity"] = within_str
                rows[idx]["cross_group_identity"] = cross_str

            families = sorted({m["family"] for m in members if m["family"]})
            group_summaries.append(
                {
                    "ortholog_id": oid,
                    "within_group_status": refined_within,
                    "cross_group_status": refined_cross,
                    "n_present": len(members),
                    "n_extracted": n_extracted,
                    "cross_group_identity": cross_str,
                    "within_group_identity": within_str,
                    "families": ";".join(families),
                    "families_consistent": str(families_consistent),
                }
            )
    finally:
        for handle in genomes.values():
            handle.close()

    out_fieldnames = [
        f
        for f in fieldnames
        if f
        not in (
            "cross_group_identity",
            "within_group_identity",
            "refined_sharing_status",
            "sharing_status",
        )
    ]
    if "within_group_status" not in out_fieldnames:
        out_fieldnames.append("within_group_status")
    if "cross_group_status" not in out_fieldnames:
        out_fieldnames.append("cross_group_status")
    out_fieldnames += ["within_group_identity", "cross_group_identity"]
    return rows, out_fieldnames, group_summaries


def is_known_family(value) -> bool:
    if value in ("", None):
        return False
    try:
        return int(value) > 0
    except (ValueError, TypeError):
        return False


def count_mixed_family_groups(rows: list[dict]) -> int:
    ortholog_groups = group_rows_by_ortholog(rows)
    n_mixed = 0
    for _, row_indices in ortholog_groups.items():
        families = {
            int(rows[idx]["family"])
            for idx in row_indices
            if rows[idx]["presence"] == "1" and is_known_family(rows[idx].get("family"))
        }
        if len(families) > 1:
            n_mixed += 1
    return n_mixed


def build_independent_insertion_events(
    rows: list[dict],
    locus_group_members_path: str | None,
) -> list[dict]:
    ortholog_groups = group_rows_by_ortholog(rows)
    events: list[dict] = []
    seen: set[tuple] = set()

    for oid, row_indices in ortholog_groups.items():
        first = rows[row_indices[0]]
        reason = first.get("cross_group_reason", "")
        cross_status = first.get("cross_group_status", "")
        within_status = first.get("within_group_status", "")

        carriers = [
            rows[idx]
            for idx in row_indices
            if rows[idx]["presence"] == "1" and is_known_family(rows[idx].get("family"))
        ]
        if not carriers:
            continue

        families = sorted({int(row["family"]) for row in carriers})
        fam_counts = Counter(int(row["family"]) for row in carriers)
        samples = sorted({row["sample"] for row in carriers})

        event_type = ""
        if reason in {"compatible_diff_family", "exact_diff_family"}:
            event_type = "same_locus_different_family"
        elif reason == "different_position":
            event_type = "different_position_insertion"
        elif len(families) > 1:
            event_type = "mixed_family_carriers"

        if not event_type:
            continue

        key = (oid, event_type, tuple(families))
        if key in seen:
            continue
        seen.add(key)

        events.append(
            {
                "ortholog_id": oid,
                "original_ortholog_id": first.get("original_ortholog_id", oid),
                "event_type": event_type,
                "cross_group_status": cross_status,
                "cross_group_reason": reason,
                "within_group_status": within_status,
                "n_carrier_families": len(families),
                "carrier_families": ",".join(str(f) for f in families),
                "carrier_family_counts": ";".join(
                    f"{fam}:{fam_counts[fam]}" for fam in families
                ),
                "present_samples": ",".join(samples),
                "locus_group_id": "",
            }
        )

    if locus_group_members_path and Path(locus_group_members_path).exists():
        import pandas as pd

        locus_df = pd.read_csv(locus_group_members_path, sep="\t")
        mixed_locus = []
        carriers = locus_df[
            (locus_df["presence"] == 1) & locus_df["family"].map(is_known_family)
        ].copy()
        if not carriers.empty:
            carriers["family"] = carriers["family"].astype(int)
            for locus_id, group in carriers.groupby("locus_group_id"):
                fam_counts = Counter(group["family"])
                if len(fam_counts) <= 1:
                    continue
                ortholog_ids = sorted(group["ortholog_id"].unique())
                mixed_locus.append(
                    {
                        "locus_group_id": locus_id,
                        "ortholog_ids": ",".join(ortholog_ids),
                        "carrier_family_counts": ";".join(
                            f"{fam}:{count}" for fam, count in sorted(fam_counts.items())
                        ),
                    }
                )

        locus_by_ortholog = {}
        for item in mixed_locus:
            for oid in item["ortholog_ids"].split(","):
                locus_by_ortholog.setdefault(oid, item["locus_group_id"])

        for event in events:
            event["locus_group_id"] = locus_by_ortholog.get(event["ortholog_id"], "")

        for item in mixed_locus:
            if any(item["locus_group_id"] == e.get("locus_group_id") for e in events if e.get("locus_group_id")):
                continue
            events.append(
                {
                    "ortholog_id": item["ortholog_ids"],
                    "original_ortholog_id": "",
                    "event_type": "locus_group_mixed_family",
                    "cross_group_status": "",
                    "cross_group_reason": "",
                    "within_group_status": "",
                    "n_carrier_families": len(item["carrier_family_counts"].split(";")),
                    "carrier_families": "",
                    "carrier_family_counts": item["carrier_family_counts"],
                    "present_samples": "",
                    "locus_group_id": item["locus_group_id"],
                }
            )

    return sorted(events, key=lambda row: (row.get("locus_group_id", ""), row["ortholog_id"]))


def write_tsv(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def resolve_until_stable(
    rows: list[dict],
    fieldnames: list[str],
    fps: dict,
    fps_by_coord: dict,
    fps_by_sample_gene: dict,
    beds_by_coord: dict,
    tandems: dict,
    cds_by_sample_gene: dict,
    codon_tolerance: int,
    max_iterations: int,
) -> tuple[list[dict], list[str], list[dict], list[dict]]:
    all_mappings: list[dict] = []
    iteration_log: list[dict] = []
    iteration = 0

    while iteration < max_iterations:
        iteration += 1
        rows, fieldnames, _ = classify_matrix(rows, fieldnames, fps, codon_tolerance)

        rows, fieldnames, within_maps, n_within = split_matrix(
            rows,
            fieldnames,
            fps_by_coord,
            fps_by_sample_gene,
            beds_by_coord,
            tandems,
            cds_by_sample_gene,
            codon_tolerance,
            cross_group_only=False,
            iteration=iteration,
        )
        all_mappings.extend(within_maps)
        iteration_log.append(
            {
                "iteration": iteration,
                "step": "within_group_split",
                "n_groups_split": n_within,
                "n_ortholog_groups": len(group_rows_by_ortholog(rows)),
                "n_rows": len(rows),
            }
        )
        if n_within > 0:
            continue

        rows, fieldnames, _ = classify_matrix(rows, fieldnames, fps, codon_tolerance)
        rows, fieldnames, cross_maps, n_cross = split_matrix(
            rows,
            fieldnames,
            fps_by_coord,
            fps_by_sample_gene,
            beds_by_coord,
            tandems,
            cds_by_sample_gene,
            codon_tolerance,
            cross_group_only=True,
            iteration=iteration,
        )
        all_mappings.extend(cross_maps)
        iteration_log.append(
            {
                "iteration": iteration,
                "step": "cross_group_split",
                "n_groups_split": n_cross,
                "n_ortholog_groups": len(group_rows_by_ortholog(rows)),
                "n_rows": len(rows),
            }
        )
        if n_cross > 0:
            continue

        iteration_log.append(
            {
                "iteration": iteration,
                "step": "stable",
                "n_groups_split": 0,
                "n_ortholog_groups": len(group_rows_by_ortholog(rows)),
                "n_rows": len(rows),
            }
        )
        break

    return rows, fieldnames, all_mappings, iteration_log


def main() -> None:
    args = parse_args()

    print("Loading auxiliary data...")
    fps_by_coord, fps_by_sample_gene = splitter.load_fingerprints(args.fingerprints)
    fps = fps_by_coord
    beds_by_coord = splitter.load_beds(args.beds)
    cds_by_sample_gene = splitter.load_gtf_cds(args.gtfs)
    tandems = splitter.load_tandems(args.tandems) if args.tandems else {}

    print(f"Reading input matrix: {args.matrix}")
    rows, fieldnames = read_matrix(args.matrix)
    print(f"  {len(rows)} rows in {len(group_rows_by_ortholog(rows))} ortholog groups")

    rows, fieldnames, all_mappings, iteration_log = resolve_until_stable(
        rows,
        fieldnames,
        fps,
        fps_by_coord,
        fps_by_sample_gene,
        beds_by_coord,
        tandems,
        cds_by_sample_gene,
        args.codon_tolerance,
        args.max_iterations,
    )

    print("\nApplying sequence identity refinement...")
    rows, fieldnames, _seq_summaries = compare_matrix(
        rows,
        fieldnames,
        args.genome_dir,
        args.flanking_length,
        args.within_group_identity_threshold,
        args.cross_group_identity_threshold,
    )

    n_mixed = count_mixed_family_groups(rows)
    n_groups = len(group_rows_by_ortholog(rows))
    note_counts = Counter(record["note"] for record in all_mappings if record["note"] != "unchanged")

    resolution_summary = list(iteration_log)
    resolution_summary.append(
        {
            "iteration": "final",
            "step": "sequence_comparison",
            "n_groups_split": 0,
            "n_ortholog_groups": n_groups,
            "n_rows": len(rows),
        }
    )
    resolution_summary.append(
        {
            "iteration": "final",
            "step": "mixed_family_ortholog_groups",
            "n_groups_split": n_mixed,
            "n_ortholog_groups": n_groups,
            "n_rows": len(rows),
        }
    )
    for note, count in sorted(note_counts.items()):
        resolution_summary.append(
            {
                "iteration": "final",
                "step": f"split_note_{note}",
                "n_groups_split": count,
                "n_ortholog_groups": n_groups,
                "n_rows": len(rows),
            }
        )

    independent_events = build_independent_insertion_events(
        rows, args.locus_group_members
    )

    write_matrix(args.output, rows, fieldnames)
    write_tsv(
        args.insertion_group_map,
        all_mappings,
        [
            "iteration",
            "resolution_pass",
            "original_id",
            "new_id",
            "n_clusters",
            "n_members",
            "n_recovered",
            "has_tandem",
            "note",
        ],
    )
    write_tsv(
        args.independent_insertion_events,
        independent_events,
        [
            "locus_group_id",
            "ortholog_id",
            "original_ortholog_id",
            "event_type",
            "cross_group_status",
            "cross_group_reason",
            "within_group_status",
            "n_carrier_families",
            "carrier_families",
            "carrier_family_counts",
            "present_samples",
        ],
    )
    write_tsv(
        args.resolution_summary,
        resolution_summary,
        ["iteration", "step", "n_groups_split", "n_ortholog_groups", "n_rows"],
    )

    print(f"\nWrote resolved matrix: {args.output}")
    print(f"Wrote insertion group map: {args.insertion_group_map}")
    print(f"Wrote independent insertion events: {args.independent_insertion_events}")
    print(f"Wrote resolution summary: {args.resolution_summary}")
    print(f"Final ortholog groups: {n_groups:,}")
    print(f"Mixed-family ortholog groups: {n_mixed:,}")
    print(f"Independent insertion events: {len(independent_events):,}")
    stable_entries = [entry for entry in iteration_log if entry["step"] == "stable"]
    n_iterations = (
        stable_entries[0]["iteration"] if stable_entries else iteration_log[-1]["iteration"]
    )
    print(f"Resolver iterations to stable: {n_iterations}")

    if n_mixed > 0:
        print("WARNING: resolver finished with mixed-family ortholog groups", file=sys.stderr)


if __name__ == "__main__":
    main()
