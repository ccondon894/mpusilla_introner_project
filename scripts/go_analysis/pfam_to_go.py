#!/usr/bin/env python3
import argparse
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests


def read_interpro2go(interpro2go):
    interpro2go_dict = defaultdict(set)
    go_pattern = re.compile(r"^GO:\d{7}$")

    with open(interpro2go) as handle:
        for line in handle:
            if line.startswith("!"):
                continue

            fields = line.strip().split()
            if not fields or not fields[0].startswith("InterPro:"):
                continue

            ipr_id = fields[0].split(":", 1)[1]
            for item in fields[1:]:
                if go_pattern.match(item):
                    interpro2go_dict[ipr_id].add(item)

    return dict(interpro2go_dict)


def collect_pfam_ids(annotation_tsv):
    df = pd.read_csv(annotation_tsv, sep="\t")
    pfam_ids = set()
    for value in df["Pfam"].dropna().astype(str):
        for pfam_id in value.split(","):
            pfam_id = pfam_id.strip()
            if pfam_id:
                pfam_ids.add(pfam_id)
    return sorted(pfam_ids)


def invert_interpro_to_pfam(interpro_to_pfam):
    pfam_to_interpro = defaultdict(set)
    for ipr_id, pfam_ids in interpro_to_pfam.items():
        if isinstance(pfam_ids, str):
            pfam_ids = [pfam_ids]
        for pfam_id in pfam_ids:
            if pfam_id:
                pfam_to_interpro[str(pfam_id).strip()].add(ipr_id)
    return dict(pfam_to_interpro)


def load_pfam_to_interpro(path):
    with open(path) as handle:
        data = json.load(handle)
    return invert_interpro_to_pfam(data)


def get_interpro_from_pfam(session, pfam_id, timeout, max_retries, retry_sleep):
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}"
    last_error = None

    for attempt in range(max_retries + 1):
        try:
            response = session.get(url, headers={"Accept": "application/json"}, timeout=timeout)
            if response.status_code == 200:
                data = response.json()
                ipr_id = data.get("metadata", {}).get("integrated")
                return {ipr_id} if ipr_id else set()
            last_error = RuntimeError(f"HTTP {response.status_code}")
        except requests.RequestException as exc:
            last_error = exc

        if attempt < max_retries:
            time.sleep(retry_sleep * (attempt + 1))

    print(f"WARNING: failed to query InterPro for {pfam_id}: {last_error}", file=sys.stderr)
    return set()


def map_pfams_to_go(pfam_ids, interpro2go, pfam_to_interpro):
    pfam_to_go = {}
    unmapped = []

    for pfam_id in pfam_ids:
        go_terms = set()
        for ipr_id in pfam_to_interpro.get(pfam_id, set()):
            go_terms.update(interpro2go.get(ipr_id, set()))
        if go_terms:
            pfam_to_go[pfam_id] = sorted(go_terms)
        else:
            pfam_to_go[pfam_id] = []
            unmapped.append(pfam_id)

    return pfam_to_go, unmapped


def main():
    parser = argparse.ArgumentParser(description="Map Pfam IDs to GO terms via InterPro.")
    parser.add_argument("interpro2go")
    parser.add_argument("annotation_tsv")
    parser.add_argument("outfile")
    parser.add_argument("--interpro-to-pfam",
                        help="Local JSON mapping InterPro IDs to Pfam IDs. Avoids InterPro API calls.")
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument("--retry-sleep", type=float, default=1.0)
    args = parser.parse_args()

    interpro2go = read_interpro2go(args.interpro2go)
    pfam_ids = collect_pfam_ids(args.annotation_tsv)
    print(f"Found {len(pfam_ids)} unique Pfam IDs", file=sys.stderr)

    pfam_to_interpro = {}
    if args.interpro_to_pfam and Path(args.interpro_to_pfam).exists():
        print(f"Loading local InterPro-to-Pfam map: {args.interpro_to_pfam}", file=sys.stderr)
        pfam_to_interpro = load_pfam_to_interpro(args.interpro_to_pfam)
    else:
        print("No local InterPro-to-Pfam map provided; querying InterPro for unique Pfams", file=sys.stderr)
        with requests.Session() as session:
            for i, pfam_id in enumerate(pfam_ids, start=1):
                print(f"Querying InterPro {i}/{len(pfam_ids)}: {pfam_id}", file=sys.stderr)
                ipr_ids = get_interpro_from_pfam(
                    session,
                    pfam_id,
                    timeout=args.timeout,
                    max_retries=args.max_retries,
                    retry_sleep=args.retry_sleep,
                )
                if ipr_ids:
                    pfam_to_interpro[pfam_id] = ipr_ids

    pfam_to_go, unmapped = map_pfams_to_go(pfam_ids, interpro2go, pfam_to_interpro)
    print(
        f"Mapped {sum(1 for terms in pfam_to_go.values() if terms)}/{len(pfam_ids)} Pfams to GO terms",
        file=sys.stderr,
    )
    if unmapped:
        print(f"Pfams without GO terms: {len(unmapped)}", file=sys.stderr)

    with open(args.outfile, "w") as handle:
        json.dump(pfam_to_go, handle, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
