#!/usr/bin/env python3
"""Fetch KEGG REST info (gene symbol, name, pathways, modules) for every
unique KO referenced by EITHER KofamScan or eggNOG-mapper.

bin/kegg_annotate.py only queries KOs that came from KofamScan hits, so any
KO eggNOG identified independently (a large fraction — eggNOG's diamond
search finds many KOs KofamScan's HMM profiles miss) never gets a symbol.
This produces a single KO-keyed lookup table covering the union of both,
reusing kegg_annotate.py's fetch/parse logic.

Output is keyed by KO (one row per unique KO), not by gene — apply it to
every KO in a gene's kegg_ko set, from whichever source, in
merge_functional_annotations.py.

Usage:
    fetch_kegg_ko_info.py --kofamscan kofamscan_result.tsv \\
        --eggnog eggnog_output.emapper.annotations \\
        --out ko_info.tsv [--cache existing_kegg_annotations.tsv]
"""
import argparse
import os
import sys
import time
from urllib.request import urlopen
from urllib.error import URLError

BATCH_SIZE = 10
SLEEP = 0.2


def fetch_ko_batch(ko_list):
    url = "https://rest.kegg.jp/get/" + "+".join(ko_list)
    for attempt in range(3):
        try:
            with urlopen(url, timeout=60) as r:
                return r.read().decode()
        except (URLError, Exception):
            if attempt == 2:
                raise
            time.sleep(5)


def parse_entries(text):
    entries = {}
    cur = None
    sec = None
    for raw in text.splitlines():
        tag = raw[:12].rstrip()
        if tag == "ENTRY":
            ko_id = raw[12:].split()[0]
            cur = {"symbols": [], "name": "", "pathways": [], "modules": []}
            entries[ko_id] = cur
            sec = None
        elif tag == "SYMBOL" and cur is not None:
            cur["symbols"] = [s.strip() for s in raw[12:].strip().split(",")]
            sec = "SYMBOL"
        elif tag == "NAME" and cur is not None:
            cur["name"] = raw[12:].strip()
            sec = "NAME"
        elif tag == "PATHWAY" and cur is not None:
            parts = raw[12:].strip().split(None, 1)
            if len(parts) == 2:
                cur["pathways"].append((parts[0], parts[1]))
            sec = "PATHWAY"
        elif tag == "MODULE" and cur is not None:
            parts = raw[12:].strip().split(None, 1)
            if len(parts) == 2:
                cur["modules"].append((parts[0], parts[1]))
            sec = "MODULE"
        elif raw.startswith("///"):
            cur = None
            sec = None
        elif raw.startswith("          ") and cur is not None:
            parts = raw.strip().split(None, 1)
            if sec == "PATHWAY" and len(parts) == 2:
                cur["pathways"].append((parts[0], parts[1]))
            elif sec == "MODULE" and len(parts) == 2:
                cur["modules"].append((parts[0], parts[1]))
        else:
            if raw and not raw[0].isspace():
                sec = None
    return entries


def kos_from_kofamscan(path):
    kos = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0] == "*":
                kos.append(parts[2])
    return kos


def kos_from_eggnog(path):
    kos = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 12 or c[11] in ('', '-'):
                continue
            for tok in c[11].split(','):
                tok = tok.strip().removeprefix('ko:')
                if tok:
                    kos.append(tok)
    return kos


def kos_already_cached(path):
    """gene-keyed kegg_annotations.tsv from kegg_annotate.py, if reusing prior fetches."""
    cached = {}
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return cached   # missing / empty cache (e.g. NO_FILE sentinel) — fetch all
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            c = line.rstrip('\n').split('\t')
            if len(c) < 8:
                continue
            ko = c[1]
            if ko and ko not in cached:
                cached[ko] = {
                    "symbols": c[2].split(',') if c[2] else [],
                    "name": c[3],
                    "pathways": list(zip(c[4].split(';'), c[5].split(';'))) if c[4] else [],
                    "modules": list(zip(c[6].split(';'), c[7].split(';'))) if c[6] else [],
                }
    return cached


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--kofamscan')
    p.add_argument('--eggnog')
    p.add_argument('--cache', help='existing gene-keyed kegg_annotations.tsv — reuse instead of re-fetching those KOs')
    p.add_argument('--out', required=True)
    args = p.parse_args()

    all_kos = []
    if args.kofamscan:
        all_kos.extend(kos_from_kofamscan(args.kofamscan))
    if args.eggnog:
        all_kos.extend(kos_from_eggnog(args.eggnog))

    seen = set()
    unique_kos = []
    for ko in all_kos:
        if ko not in seen:
            seen.add(ko)
            unique_kos.append(ko)
    sys.stderr.write(f"Unique KOs referenced (KofamScan + eggNOG): {len(unique_kos)}\n")

    ko_info = {}
    to_fetch = unique_kos
    if args.cache:
        ko_info = kos_already_cached(args.cache)
        to_fetch = [k for k in unique_kos if k not in ko_info]
        sys.stderr.write(f"Already cached: {len(ko_info)}  -  still need to fetch: {len(to_fetch)}\n")

    sys.stderr.write(f"Fetching {len(to_fetch)} KOs in batches of {BATCH_SIZE}...\n")
    for i in range(0, len(to_fetch), BATCH_SIZE):
        batch = to_fetch[i:i + BATCH_SIZE]
        text = fetch_ko_batch(batch)
        ko_info.update(parse_entries(text))
        sys.stderr.write(f"  {min(i + BATCH_SIZE, len(to_fetch))}/{len(to_fetch)}\r")
        time.sleep(SLEEP)
    sys.stderr.write("\n")

    header = "\t".join(["ko", "gene_symbols", "ko_name", "pathway_ids", "pathway_names", "module_ids", "module_names"])
    with open(args.out, "w") as out:
        out.write(header + "\n")
        for ko in unique_kos:
            info = ko_info.get(ko, {})
            symbols = ",".join(info.get("symbols", []))
            name = info.get("name", "")
            pathway_ids = ";".join(p[0] for p in info.get("pathways", []))
            pathway_names = ";".join(p[1] for p in info.get("pathways", []))
            module_ids = ";".join(m[0] for m in info.get("modules", []))
            module_names = ";".join(m[1] for m in info.get("modules", []))
            out.write("\t".join([ko, symbols, name, pathway_ids, pathway_names, module_ids, module_names]) + "\n")

    sys.stderr.write(f"Written: {args.out}\n")


if __name__ == "__main__":
    main()
