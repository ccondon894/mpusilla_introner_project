#!/usr/bin/env python3
import argparse
import json
import sys
import time

import pandas as pd
import requests


KEGG_GET_LIMIT = 10


def collect_ko_ids(annotation_tsv):
    df = pd.read_csv(annotation_tsv, sep="\t")
    ko_ids = set()
    for value in df["KO"].dropna().astype(str):
        for ko_id in value.split(","):
            ko_id = ko_id.strip()
            if ko_id:
                ko_ids.add(ko_id)
    return sorted(ko_ids)


def chunks(items, size):
    for start in range(0, len(items), size):
        yield items[start:start + size]


def process_line_for_go_terms(line, go_terms):
    def normalize_go(raw_go):
        raw_go = raw_go.strip()
        if not raw_go:
            return None
        return raw_go if raw_go.startswith("GO:") else f"GO:{raw_go}"

    parts = line.strip().split()
    for i, part in enumerate(parts):
        if part == "GO:":
            go_terms.update(
                go_id for go_id in (normalize_go(go) for go in parts[i + 1:])
                if go_id
            )
        elif part == "GO" and i + 1 < len(parts) and parts[i + 1] == ":":
            go_terms.update(
                go_id for go_id in (normalize_go(go) for go in parts[i + 2:])
                if go_id
            )
        elif part.startswith("GO:"):
            go_terms.add(part)


def parse_kegg_entries(text):
    ko_to_go = {}
    for entry in text.split("///"):
        entry = entry.strip()
        if not entry:
            continue

        ko_id = None
        go_terms = set()
        in_dblinks_section = False

        for line in entry.splitlines():
            if line.startswith("ENTRY"):
                parts = line.split()
                if len(parts) >= 2:
                    ko_id = parts[1]
            elif line.startswith("DBLINKS"):
                in_dblinks_section = True
                process_line_for_go_terms(line, go_terms)
            elif in_dblinks_section and line.startswith("            "):
                process_line_for_go_terms(line, go_terms)
            elif in_dblinks_section:
                in_dblinks_section = False

        if ko_id:
            ko_to_go[ko_id] = go_terms

    return ko_to_go


def fetch_batch(session, ko_ids, timeout, max_retries, retry_sleep):
    url = "https://rest.kegg.jp/get/" + "+".join(ko_ids)
    last_error = None

    for attempt in range(max_retries + 1):
        try:
            response = session.get(url, timeout=timeout)
            if response.status_code == 200:
                return parse_kegg_entries(response.text)
            last_error = RuntimeError(f"HTTP {response.status_code}")
        except requests.RequestException as exc:
            last_error = exc

        if attempt < max_retries:
            time.sleep(retry_sleep * (attempt + 1))

    print(
        f"WARNING: failed to query KEGG batch {','.join(ko_ids)}: {last_error}",
        file=sys.stderr,
    )
    return {}


def main():
    parser = argparse.ArgumentParser(description="Map KEGG Orthology IDs to GO terms.")
    parser.add_argument("annotation_tsv")
    parser.add_argument("outfile")
    parser.add_argument("--batch-size", type=int, default=KEGG_GET_LIMIT,
                        help=f"Number of KO IDs per KEGG get request, max {KEGG_GET_LIMIT}.")
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument("--retry-sleep", type=float, default=1.0)
    args = parser.parse_args()

    batch_size = min(max(1, args.batch_size), KEGG_GET_LIMIT)
    ko_ids = collect_ko_ids(args.annotation_tsv)
    ko_dict = {ko_id: set() for ko_id in ko_ids}

    print(f"Found {len(ko_ids)} unique KO IDs", file=sys.stderr)
    with requests.Session() as session:
        for i, batch in enumerate(chunks(ko_ids, batch_size), start=1):
            print(f"Querying KEGG batch {i}: {len(batch)} KO IDs", file=sys.stderr)
            batch_results = fetch_batch(
                session,
                batch,
                timeout=args.timeout,
                max_retries=args.max_retries,
                retry_sleep=args.retry_sleep,
            )
            for ko_id, go_terms in batch_results.items():
                if ko_id in ko_dict:
                    ko_dict[ko_id].update(go_terms)

    ko_dict_for_json = {
        ko_id: sorted(go_terms)
        for ko_id, go_terms in ko_dict.items()
    }

    with open(args.outfile, "w") as handle:
        json.dump(ko_dict_for_json, handle, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
