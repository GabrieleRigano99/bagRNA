#!/usr/bin/env python3
"""
Fetch KEGG module definitions for every module found in this genome's
kegg_annotations.tsv (from bin/kegg_annotate.py), and score completeness
per pathway STEP rather than per raw KO: each line of a module's ORTHOLOGY
block is one step (a reaction), and a step counts as satisfied if the
genome has ANY ONE of that step's alternative KOs or EC numbers — matching
how these modules are actually meant to be read (KEGG lists 2-6 isozyme
KOs per step for well-studied pathways; a real organism normally carries
only one of them, so scoring by raw KO count structurally penalizes
exactly the universal, heavily-annotated pathways like glycolysis).

Evidence comes from two sources: KofamScan/DeepKOALA's kegg_annotations.tsv
(KOs) and, when given, eggNOG's own KEGG_ko and EC columns — eggNOG often
calls genes the KO-classifier pipeline misses, so skipping it understates
real completeness for genes that only eggNOG annotated.

Usage: fetch_kegg_module_completeness.py kegg_annotations.tsv output.json [--eggnog eggnog.emapper.annotations]
"""
import argparse
import json
import re
import sys
import time
from urllib.request import urlopen

BATCH_SIZE = 10
SLEEP = 0.3
KO_RE = re.compile(r"K\d{5}")
EC_RE = re.compile(r"\d+\.\d+\.\d+\.(?:\d+|-)")


def collect_found(kegg_annotations_path, eggnog_path):
    """Modules relevant to this genome are discovered from kegg_annotations.tsv
    alone (unchanged scope); found_kos/found_ecs — the evidence used to score
    each step — are the union of that file and eggNOG's independent calls."""
    kos_by_mod = {}
    mod_names = {}
    found_kos = set()
    with open(kegg_annotations_path) as f:
        next(f, None)
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            ko = cols[1]
            found_kos.add(ko)
            mod_ids = cols[6]
            if not mod_ids:
                continue
            ids = mod_ids.split(";")
            names = cols[7].split(";") if cols[7] else []
            for i, m in enumerate(ids):
                kos_by_mod.setdefault(m, set())
                if i < len(names):
                    mod_names[m] = names[i]

    found_ecs = set()
    if eggnog_path:
        with open(eggnog_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12:
                    continue
                ec_field, ko_field = cols[10], cols[11]
                if ec_field and ec_field != "-":
                    found_ecs.update(e.replace("ec:", "") for e in ec_field.split(","))
                if ko_field and ko_field != "-":
                    found_kos.update(k.replace("ko:", "") for k in ko_field.split(","))

    return sorted(kos_by_mod), mod_names, found_kos, found_ecs


def fetch_batch(ids):
    url = "https://rest.kegg.jp/get/" + "+".join(ids)
    for attempt in range(3):
        try:
            with urlopen(url, timeout=60) as r:
                return r.read().decode()
        except Exception:
            if attempt == 2:
                raise
            time.sleep(5)


def parse_modules(text):
    """Each ORTHOLOGY line becomes one step: {kos, ecs} sets of alternatives.
    Scoping (which section a continuation line belongs to) matches
    kegg_annotate.py's KO-entry parser — a bare column-prefix match would
    also pick up unrelated tokens from other sections (e.g. COMPLETE)."""
    entries = {}
    cur_id = None
    sec = None
    steps = []

    def _step_from(raw):
        kos, ecs = KO_RE.findall(raw), EC_RE.findall(raw)
        return {"kos": set(kos), "ecs": set(ecs)} if (kos or ecs) else None

    for raw in text.splitlines():
        tag = raw[:12].rstrip()
        if tag == "ENTRY":
            if cur_id is not None:
                entries[cur_id] = steps
            cur_id = raw[12:].split()[0]
            steps = []
            sec = None
        elif raw.startswith("///"):
            if cur_id is not None:
                entries[cur_id] = steps
            cur_id = None
            steps = []
            sec = None
        elif tag:
            sec = tag
            if sec == "ORTHOLOGY":
                step = _step_from(raw)
                if step:
                    steps.append(step)
        elif raw.startswith("            ") and sec == "ORTHOLOGY":
            step = _step_from(raw)
            if step:
                steps.append(step)
    if cur_id is not None and cur_id not in entries:
        entries[cur_id] = steps
    return entries


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("kegg_annotations")
    p.add_argument("output_json")
    p.add_argument("--eggnog", default=None, help="eggNOG emapper.annotations (optional, adds KO/EC evidence)")
    args = p.parse_args()

    module_ids, mod_names, found_kos, found_ecs = collect_found(args.kegg_annotations, args.eggnog)
    print(f"{len(module_ids)} unique modules found in this genome's KEGG annotation "
          f"({len(found_kos)} KOs, {len(found_ecs)} EC numbers as evidence)", file=sys.stderr)

    total_required = {}
    for i in range(0, len(module_ids), BATCH_SIZE):
        batch = module_ids[i:i + BATCH_SIZE]
        text = fetch_batch(batch)
        total_required.update(parse_modules(text))
        print(f"  {min(i + BATCH_SIZE, len(module_ids))}/{len(module_ids)}", end="\r", file=sys.stderr)
        time.sleep(SLEEP)
    print(file=sys.stderr)

    results = []
    for m in module_ids:
        steps = total_required.get(m)
        if not steps:
            continue
        satisfied = sum(1 for s in steps if (s["kos"] & found_kos) or (s["ecs"] & found_ecs))
        pct = 100.0 * satisfied / len(steps)
        results.append({
            "module_id": m,
            "module_name": mod_names.get(m, m),
            "found": satisfied,
            "total": len(steps),
            "pct": round(pct, 1),
        })
    results.sort(key=lambda x: (-x["pct"], -x["total"]))

    with open(args.output_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Wrote {args.output_json} ({len(results)} modules)", file=sys.stderr)
    print("\nFully complete (100%) modules:", file=sys.stderr)
    for r in results:
        if r["pct"] == 100.0:
            print(f"  {r['module_id']}  {r['module_name']}  ({r['found']}/{r['total']} steps)", file=sys.stderr)


if __name__ == "__main__":
    main()
