#!/usr/bin/env python3
"""
Enrich KofamScan TSV with KEGG gene symbols, pathways, and modules.
Queries rest.kegg.jp/get in batches of 10 KOs.

Input : kofamscan_result.tsv  (tab-sep: * gene_id KO threshold score evalue description)
Output: kegg_annotations.tsv  (tab-sep: gene_id KO gene_symbols ko_name
                                         pathway_ids pathway_names module_ids module_names)
"""
import sys
import time
from urllib.request import urlopen
from urllib.error import URLError

BATCH_SIZE = 10
SLEEP      = 0.2   # seconds between requests


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
    """
    Parse KEGG flat-file text into:
        { KO_id: { symbols, name, pathways: [(id, name)], modules: [(id, name)] } }
    """
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
            # continuation line — only relevant for PATHWAY and MODULE
            parts = raw.strip().split(None, 1)
            if sec == "PATHWAY" and len(parts) == 2:
                cur["pathways"].append((parts[0], parts[1]))
            elif sec == "MODULE" and len(parts) == 2:
                cur["modules"].append((parts[0], parts[1]))

        else:
            if raw and not raw[0].isspace():
                sec = None

    return entries


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: kegg_annotate.py kofamscan_result.tsv output.tsv")

    kofam_path, out_path = sys.argv[1], sys.argv[2]

    # Read kofamscan: col0=*, col1=gene_id, col2=KO
    rows = []
    with open(kofam_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0] == "*":
                rows.append((parts[1], parts[2]))   # gene_id, KO

    # Unique KOs preserving order
    seen = set()
    unique_kos = []
    for _, ko in rows:
        if ko not in seen:
            seen.add(ko)
            unique_kos.append(ko)

    print(f"Fetching {len(unique_kos)} unique KOs in batches of {BATCH_SIZE}...",
          file=sys.stderr)

    ko_info = {}
    for i in range(0, len(unique_kos), BATCH_SIZE):
        batch = unique_kos[i:i + BATCH_SIZE]
        text  = fetch_ko_batch(batch)
        ko_info.update(parse_entries(text))
        print(f"  {min(i + BATCH_SIZE, len(unique_kos))}/{len(unique_kos)}",
              file=sys.stderr, end="\r")
        time.sleep(SLEEP)

    print(file=sys.stderr)

    header = "\t".join([
        "gene_id", "ko",
        "gene_symbols", "ko_name",
        "pathway_ids", "pathway_names",
        "module_ids",  "module_names",
    ])

    with open(out_path, "w") as out:
        out.write(header + "\n")
        for gene_id, ko in rows:
            info = ko_info.get(ko, {})
            symbols      = ",".join(info.get("symbols", []))
            name         = info.get("name", "")
            pathway_ids  = ";".join(p[0] for p in info.get("pathways", []))
            pathway_names = ";".join(p[1] for p in info.get("pathways", []))
            module_ids   = ";".join(m[0] for m in info.get("modules", []))
            module_names = ";".join(m[1] for m in info.get("modules", []))
            out.write("\t".join([
                gene_id, ko,
                symbols, name,
                pathway_ids, pathway_names,
                module_ids, module_names,
            ]) + "\n")

    print(f"Written: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
