"""GTF interval-based gene annotation helpers for RF feature matrices."""

from __future__ import annotations

import bisect
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


GENE_FEATURE_PRIORITY = ("gene", "transcript", "mRNA", "exon", "CDS")


def parse_gtf_attrs(attr_text: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]+)"', attr_text))


@dataclass(frozen=True)
class GeneHit:
    gene: str | pd.NA
    source: str | pd.NA
    overlap_bp: int
    n_overlaps: int


class GtfGeneAnnotator:
    """Annotate genomic intervals with gene IDs from a GTF."""

    def __init__(self, gtf_path, feature_priority=GENE_FEATURE_PRIORITY):
        self.gtf_path = Path(gtf_path)
        self.feature_priority = tuple(feature_priority)
        self.intervals = self._read_intervals()

    def _read_intervals(self):
        intervals = {
            feature_type: defaultdict(list)
            for feature_type in self.feature_priority
        }
        with open(self.gtf_path) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue
                contig, _, feature_type, start, end, _, _, _, attr_text = fields
                if feature_type not in intervals:
                    continue
                attrs = parse_gtf_attrs(attr_text)
                gene = attrs.get("gene_id")
                if not gene:
                    continue
                # GTF coordinates are 1-based closed; use 0-based half-open.
                intervals[feature_type][contig].append((int(start) - 1, int(end), gene))

        for feature_type in intervals:
            for contig, rows in intervals[feature_type].items():
                rows.sort(key=lambda row: (row[0], row[1], row[2]))
                starts = [row[0] for row in rows]
                intervals[feature_type][contig] = {"starts": starts, "rows": rows}
        return intervals

    def best_hit(self, contig, start, end) -> GeneHit:
        if pd.isna(contig) or pd.isna(start) or pd.isna(end):
            return GeneHit(pd.NA, pd.NA, 0, 0)

        start = int(start)
        end = int(end)
        if end <= start:
            end = start + 1

        for feature_type in self.feature_priority:
            contig_index = self.intervals[feature_type].get(str(contig))
            if not contig_index:
                continue
            hits = self._overlaps(contig_index, start, end)
            if not hits:
                continue
            ranked = sorted(
                hits,
                key=lambda hit: (
                    -hit["overlap_bp"],
                    hit["interval_start"],
                    hit["interval_end"],
                    hit["gene"],
                ),
            )
            best = ranked[0]
            return GeneHit(
                best["gene"],
                feature_type,
                int(best["overlap_bp"]),
                len({hit["gene"] for hit in hits}),
            )
        return GeneHit(pd.NA, pd.NA, 0, 0)

    @staticmethod
    def _overlaps(contig_index, start, end):
        starts = contig_index["starts"]
        rows = contig_index["rows"]
        idx = bisect.bisect_right(starts, start)
        hits = []

        scan_idx = max(0, idx - 1)
        while scan_idx > 0 and rows[scan_idx - 1][1] > start:
            scan_idx -= 1

        while scan_idx < len(rows):
            interval_start, interval_end, gene = rows[scan_idx]
            if interval_start >= end:
                break
            overlap = min(end, interval_end) - max(start, interval_start)
            if overlap > 0:
                hits.append({
                    "gene": gene,
                    "interval_start": interval_start,
                    "interval_end": interval_end,
                    "overlap_bp": overlap,
                })
            scan_idx += 1
        return hits

    def annotate_frame(
        self,
        df: pd.DataFrame,
        contig_col="contig",
        start_col="feature_start",
        end_col="feature_end",
        gene_col="gene",
    ) -> pd.DataFrame:
        out = df.copy()
        if gene_col in out.columns and "gene_from_genotype_matrix" not in out.columns:
            out["gene_from_genotype_matrix"] = out[gene_col]

        hits = [
            self.best_hit(contig, start, end)
            for contig, start, end in zip(
                out[contig_col],
                out[start_col],
                out[end_col],
            )
        ]
        out[gene_col] = [hit.gene for hit in hits]
        out["gtf_gene_source"] = [hit.source for hit in hits]
        out["gtf_gene_overlap_bp"] = [hit.overlap_bp for hit in hits]
        out["gtf_gene_n_overlaps"] = [hit.n_overlaps for hit in hits]
        return out
