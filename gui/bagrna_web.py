#!/usr/bin/env python3
"""Browser GUI for the bagRNA pipeline, served on localhost. Standard library
only. Usage:  python3 gui/bagrna_web.py [--port 8765] [--no-browser]
Binds 127.0.0.1 only: the server launches processes, never expose it.
"""
import argparse
import json
import os
import re
import shlex
import signal
import subprocess
import threading
import time
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs

PIPELINE_DIR = Path(__file__).resolve().parent.parent
HTML = Path(__file__).with_name("bagrna_web.html")
# Launch environment required by bagRNA (see CLAUDE.md "Running the pipeline").
NXF_VER = "25.04.6"
JAVA_HOME = "/usr/lib/jvm/java-21-openjdk-amd64"


# ---- form fields: nextflow_schema.json (labels/groups/enums) + nextflow.config (truth) ----
# The schema lags the pipeline (no IPS6_databases_path, eggnog_data_dir, lr_manifest,
# SETUP flags...; wrong defaults for scoring/orientation), so every `params {}` entry in
# nextflow.config becomes a field and the schema only decorates it.
MAIN_GROUPS = ("Required", "Databases & tools", "Output options")
DB_FIELDS = {"databases", "eggnog_data_dir", "IPS6_databases_path", "signalp_path", "phobius_path",
             "dbcan_db", "merops_db", "phi_base_db"}
BROWSE_BY_NAME = {"eggnog_data_dir": "dir", "IPS6_databases_path": "dir", "dbcan_db": "dir",
                  "db_dir": "dir", "outdir": "dir", "star_manifest": "files", "pfam_db": "file",
                  "lr_manifest": "file", "merops_db": "file", "phi_base_db": "file", "go_obo": "file",
                  "gene2product": "file", "id_map": "file", "scoring": "file"}
MULTI = {"star_manifest"}  # comma-separated in the form -> repeated flag ("Repeat flag for multiple manifests")


def load_fields():
    schema = json.loads((PIPELINE_DIR / "nextflow_schema.json").read_text())
    sprops = {}
    for gname, g in (schema.get("definitions") or schema.get("$defs", {})).items():
        for pname, p in g.get("properties", {}).items():
            if not p.get("hidden"):
                sprops[pname] = (p, "Required" if gname == "mandatory_inputs" else g.get("title", gname))
    cfg = (PIPELINE_DIR / "nextflow.config").read_text()
    block = re.search(r"^params\s*\{(.*?)^\}", cfg, re.S | re.M).group(1)
    fields, browse, section = [], {}, ""
    for line in block.splitlines():
        h = re.match(r"\s*//\s*─+\s*(.*?)\s*─+\s*$", line)
        if h:
            section = h.group(1)
            continue
        m = re.match(r"\s*(\w+)\s*=\s*(\"[^\"]*\"|'[^']*'|[^\s/]+)\s*(?://\s*(.*))?$", line)
        if not m or m.group(1) in ("help", "tracedir"):
            continue
        name, raw, comment = m.group(1), m.group(2), m.group(3) or ""
        p, sgroup = sprops.get(name, ({}, None))
        if raw in ("true", "false"):
            default = raw == "true"
        else:
            default = "" if raw == "null" else raw.strip("\"'").replace("${projectDir}", str(PIPELINE_DIR))
        if p.get("type") == "boolean" and not isinstance(default, bool):
            default = False
        kind = "bool" if isinstance(default, bool) else (p["enum"] if "enum" in p else "text")
        if isinstance(kind, list) and default not in kind:
            kind = [""] + kind
        group = "Databases & tools" if name in DB_FIELDS else sgroup or section
        if group not in MAIN_GROUPS:
            group = "Advanced: " + group
        mode = "setup" if "-entry SETUP" in section else "run"
        fields.append((name, name, default, kind, group, mode, p.get("description") or comment or name))
        b = BROWSE_BY_NAME.get(name) or {"file-path": "file", "directory-path": "dir"}.get(p.get("format"))
        if b:
            browse[name] = b
    profiles = re.findall(r"^    (\w+)\s*\{", re.search(r"^profiles\s*\{(.*?)^\}", cfg, re.S | re.M).group(1), re.M)
    return fields, browse, profiles


FIELDS, BROWSE, PROFILES = load_fields()


def list_dir(path):
    p = Path(path or Path.home()).expanduser()
    if p.is_file():
        p = p.parent
    while not p.is_dir() and p != p.parent:
        p = p.parent
    out = []
    try:
        for e in sorted(p.iterdir(), key=lambda e: (not e.is_dir(), e.name.lower())):
            if not e.name.startswith("."):
                out.append({"name": e.name, "dir": e.is_dir()})
    except OSError:
        pass
    return {"path": str(p.resolve()), "parent": str(p.resolve().parent), "entries": out}


def build_command(p):
    """p: dict of field values + entry/profile/resume/stub. Returns (args, errors)."""
    g = lambda k: str(p.get(k) or "").strip()
    setup = p.get("entry") == "setup"
    errors = []
    # Same rules as check_required() in main.nf; Nextflow stays the authority.
    if setup:
        need = ["db_dir"]
    elif p.get("functional_anno_only"):
        need = ["genome_fasta", "species", "protein_fasta", "submission_template"]
    else:
        need = ["genome_fasta", "species", "prot_evidence", "busco_lineage", "submission_template"]
        if not (p.get("no_helixer") or g("helixer_gff") or g("helixer_lineage")):
            errors.append("Helixer input required: helixer_gff, helixer_lineage, or no_helixer.")
    errors += [f"{k} is required." for k in need if not g(k)]
    if errors:
        return None, errors

    args = ["nextflow", "run", "main.nf"] + (["-entry", "SETUP"] if setup else [])
    for name, _, default, kind, _, mode, _ in FIELDS:
        if mode != ("setup" if setup else "run"):
            continue
        if kind == "bool":
            on = bool(p.get(name, default))
            if on != default:  # default-true flags (jaccard_clip...) need an explicit "false"
                args += [f"--{name}"] + ([] if on else ["false"])
        elif g(name) and g(name) != default:
            for v in (g(name).split(",") if name in MULTI else [g(name)]):
                args += [f"--{name}", os.path.expanduser(v.strip())]
    if g("profile") not in ("", "standard"):
        args += ["-profile", g("profile")]
    if p.get("resume"):
        args.append("-resume")
    if p.get("stub"):
        args.append("-stub-run")
    return args, []


# Nextflow's own progress-table line, e.g. "[2b/131140] STR…:COMPLEASM_GENOME (genome) | 1 of 1, cached: 1 ✔"
STATUS_RE = re.compile(
    r"^\[(?P<hash>[^\]]*)\]\s+(?P<name>\S+)(?:\s+\((?P<tag>[^)]*)\))?\s*\|\s*"
    r"(?P<done>\d+)\s+of\s+(?P<total>\d+)(?:,\s*(?P<extra>[\w\s:]+))?\s*(?P<check>✔)?\s*$"
)
PEND_RE = re.compile(r"^\[-\s*\]\s+(?P<name>\S+)\s+-\s*$")  # process listed, nothing submitted yet
FAILED_RE = re.compile(r"Error executing process\s*>\s*'([^']+)'")


def _sources():
    for d in ("modules", "workflows", "subworkflows"):
        yield from (PIPELINE_DIR / d).rglob("*.nf")
    yield PIPELINE_DIR / "main.nf"


# Process names from the pipeline's own sources (incl. `include { X as Y }` aliases), never hand-typed.
KNOWN = sorted({n for f in _sources() for n in
                re.findall(r"^process\s+(\w+)", f.read_text(), re.M) + re.findall(r"\bas\s+([A-Z][A-Z0-9_]+)\b", f.read_text())})


def _stages():
    out = {}
    for f, st in (("download_databases", "Setup"), ("functional_annotation", "Functional"), ("structural_annotation", "Structural")):
        txt = (PIPELINE_DIR / "workflows" / f"{f}.nf").read_text()
        for blk in re.findall(r"include\s*\{([^}]*)\}", txt):
            for item in re.split(r"[;\n]", blk):
                if item.strip():
                    out[item.split()[-1]] = st
    for f in (PIPELINE_DIR / "subworkflows").rglob("*.nf"):  # InterProScan 6's own processes
        out.update({n: "InterProScan 6" for n in re.findall(r"^process\s+(\w+)|\bas\s+(\w+)", f.read_text(), re.M) for n in n if n})
    return out


WF_OF = _stages()


def _tags():  # static `tag "literal"` of each process, to resolve names Nextflow cut to e.g. "STR…PTS (busco_protein_scores)"
    out = {}
    for f in _sources():
        for blk in re.split(r"^(?=process\s)", f.read_text(), flags=re.M):
            m = re.match(r"process\s+(\w+)", blk)
            if m:
                out[m[1]] = set(re.findall(r"^\s*tag\s+[\"']([^\"'$]+)[\"']", blk, re.M))
    return out


TAGS = _tags()


def canon(name, tag=""):
    """Nextflow ellipsizes 'WORKFLOW:SUB:PROCESS' into 'pre…post' (truncation varies between
    reprints, and can leave only 'STR…'). Return the process name (last segment); the tag
    breaks ties. Unresolvable names are returned unchanged."""
    if "(" in name:  # tag glued to a cut name: "FUN…(kegg_module_completeness)"
        name, _, t = name.partition("(")
        tag = tag or t.rstrip(")")
    if "…" not in name:
        return name.rsplit(":", 1)[-1]
    post = name.partition("…")[2]
    if post.endswith(")"):  # the cut ate the whole name and "(": what is left is the tag
        tag, post = tag or post[:-1], ""
    last = post.rsplit(":", 1)[-1]
    full = ":" in post  # the segment after the last ':' is a complete process name
    cands = [k for k in KNOWN if (k == last if full else k.endswith(last))]
    if len(cands) > 1 and tag:
        t = tag.lower()
        tagged = lambda k: k.lower() == t or tag in TAGS.get(k, ())
        cands = [k for k in cands if tagged(k)] or [k for k in cands if k.lower() in t or t in k.lower()] or cands
    return cands[0] if len(cands) == 1 else name


def read_failure(workdir):
    """Tail of the failed task's stderr and its command, from the work dir
    Nextflow printed in its own error message."""
    d = Path(workdir)
    out = {"err": "", "cmd": ""}
    for key, name, tail in (("err", ".command.err", True), ("cmd", ".command.sh", False)):
        try:
            lines = (d / name).read_text(errors="replace").splitlines()
            out[key] = "\n".join(lines[-40:] if tail else lines[:40])
        except OSError:
            pass
    return out


def _outdir(outdir):
    return (PIPELINE_DIR / Path(outdir or "results").expanduser()).resolve()  # relative = relative to pipeline dir (nextflow's cwd)


def results(outdir):
    """Parse a finished output dir. Primary source: the pipeline's own annotation_report/*_annotation_summary.tbl."""
    o = _outdir(outdir)
    rep = o / "annotation_report"
    tbls = sorted(rep.glob("*_annotation_summary.tbl")) if rep.is_dir() else []
    groups = {}
    for line in (tbls[0].read_text().splitlines()[1:] if tbls else []):
        c = line.split("\t")
        if len(c) == 3:
            groups.setdefault(c[0], []).append(c[1:])
    if not tbls:  # run stopped before the report stage: same groups, from the raw files
        comp = o / "compleasm_genome" / "compleasm_out" / "summary.txt"
        if comp.is_file():
            t = comp.read_text()
            pc = dict(re.findall(r"^([SDFIM]):([\d.]+)%", t, re.M))
            groups["genome"] = [["busco_lineage", (re.search(r"lineage: (\S+)", t) or [0, ""])[1]],
                                ["busco_complete_single_pct", pc.get("S")], ["busco_complete_duplicated_pct", pc.get("D")],
                                ["busco_fragmented_pct", pc.get("F")], ["busco_missing_pct", pc.get("M")],
                                ["busco_total_groups", (re.search(r"^N:(\d+)", t, re.M) or [0, ""])[1]]]
        gffs = sorted((o / "structural_annotation").glob("final_struct_*.gff"))
        if gffs:
            cnt = {}
            for ln in open(gffs[0]):
                f = ln.split("\t")
                if len(f) > 8 and f[2] in ("gene", "mRNA", "ncRNA", "tRNA", "rRNA"):
                    cnt[f[2]] = cnt.get(f[2], 0) + 1
            groups["structural_feature"] = [[k, str(v)] for k, v in cnt.items()]
        if not groups:
            return {"found": False, "outdir": str(o)}
    files = [f.name for f in sorted(rep.iterdir()) if f.suffix in (".html", ".pdf", ".tbl")] if rep.is_dir() else []
    return {"found": True, "outdir": str(o), "groups": groups, "files": files,
            "trace": (o / "pipeline_info" / "execution_report.html").is_file()}


def open_file(outdir, rel):
    """xdg-open a file or the folder, only inside the output dir."""
    o = _outdir(outdir)
    t = (o / rel).resolve() if rel else o
    if (t == o or o in t.parents) and t.exists():
        subprocess.Popen(["xdg-open", str(t)])
        return True
    return False


class Run:
    """Single pipeline run; state is polled by the page."""

    def __init__(self):
        self.lock = threading.Lock()
        self.reset()

    def reset(self):
        self.proc = None
        self.log = []          # list of [text, is_error]
        self.status = "idle"   # idle|running|done|failed
        self.code = None
        self.cmd = ""
        self.t0 = None
        self.outdir = ""
        self.stopped = False
        self.rows = {}         # process -> dict
        self.order = []
        self.failed = set()
        self.times = {}        # process -> [start, end|None]  (epoch seconds)
        self.failures = []     # [{"proc","workdir"}]
        self._last_fail = ""
        self._wd_next = False
        self._err_cache = {}

    def start(self, params):
        with self.lock:
            if self.status == "running":
                return ["A run is already in progress."]
            args, errors = build_command(params)
            if errors:
                return errors
            env = os.environ.copy()
            for k, key in (("NXF_VER", "nxf_ver"), ("JAVA_HOME", "java_home")):
                if str(params.get(key) or "").strip():
                    env[k] = params[key].strip()
            self.reset()
            self.cmd = shlex.join(args)
            self.log.append(["$ " + self.cmd + "\n", False])
            try:
                # own session: Stop must reach nextflow's whole process group
                self.proc = subprocess.Popen(
                    args, cwd=PIPELINE_DIR, env=env, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, text=True, bufsize=1, start_new_session=True,
                )
            except FileNotFoundError:
                self.status = "idle"
                return ["Could not launch `nextflow`. Is it on your PATH?"]
            self.status = "running"
            self.outdir = str(params.get("outdir", "")).strip()
            self.t0 = time.time()
            threading.Thread(target=self._read, args=(self.proc,), daemon=True).start()
        return []

    def _read(self, proc):
        for line in proc.stdout:
            with self.lock:
                self._line(line)
        proc.wait()
        with self.lock:
            self.code = proc.returncode
            self.status = "done" if proc.returncode == 0 else "stopped" if self.stopped else "failed"

    def _line(self, line):
        err = line.startswith("ERROR") or "Error executing process" in line
        self.log.append([line if line.endswith("\n") else line + "\n", err])
        m = FAILED_RE.search(line)
        if m:
            self._last_fail = m.group(1)
            self.failed.add(canon(*self._split(m.group(1))))
        if line.strip() == "Work dir:":
            self._wd_next = True
        elif self._wd_next and line.strip():
            self._wd_next = False
            if all(f["workdir"] != line.strip() for f in self.failures):  # nextflow prints the error twice
                self.failures.append({"proc": self._last_fail, "workdir": line.strip()})
        p = PEND_RE.match(line.strip())
        s = STATUS_RE.match(line.strip())
        if p:
            self.log[-1].append(canon(p["name"]))
            self._row(self.log[-1][2], "", "", 0, 0, "", "")
        elif s:
            d = s.groupdict()
            self.log[-1].append(canon(d["name"], d["tag"] or ""))
            self._row(self.log[-1][2], d["tag"] or "", d["hash"].strip(),
                      int(d["done"]), int(d["total"]), (d["extra"] or "").strip(), d["name"])

    @staticmethod
    def _split(full):  # "WF:PROC (tag)" -> ("WF:PROC", "tag")
        n, _, t = full.partition(" (")
        return n, t.rstrip(")")

    def _row(self, name, tag, hsh, done, total, extra, raw):
        if name not in self.rows:
            self.order.append(name)
        self.rows[name] = {"name": name, "tag": tag, "done": done, "total": total, "hash": hsh,
                           "extra": extra, "stage": WF_OF.get(name, "Other")}
        if total:
            t = self.times.setdefault(name, [time.time(), None])
            if t[1] is None and total <= done:
                t[1] = time.time()

    def stop(self):
        with self.lock:
            if self.proc and self.status == "running":
                self.stopped = True
                os.killpg(self.proc.pid, signal.SIGTERM)
                self.log.append(["\n[GUI] Sent SIGTERM to nextflow's process group.\n", True])

    def snapshot(self, since):
        with self.lock:
            rows, t0 = [], self.t0 or 0
            for n in self.order:
                r = dict(self.rows[n])
                if n in self.failed: st = "failed"
                elif r["hash"] == "skipped" or "cached" in r["extra"] or "stored" in r["extra"]: st = "cached"
                elif r["total"] and r["done"] >= r["total"]: st = "done"
                elif r["total"]: st = "running"
                else: st = "pending"
                r["state"] = st
                a, b = self.times.get(n, [None, None])
                r["start"] = round(a - t0, 1) if a else None
                r["end"] = round(b - t0, 1) if b else None
                rows.append(r)
            fails = []
            for f in self.failures:
                if f["workdir"] not in self._err_cache:
                    self._err_cache[f["workdir"]] = read_failure(f["workdir"])
                fl = {**f, **self._err_cache[f["workdir"]]}
                if not fl["err"]:  # no .command.err (native/exec process, missing outputs): show nextflow's own message
                    txt = [l[0].rstrip() for l in self.log]
                    i = next((i for i, l in enumerate(txt) if "Error executing process" in l and f["proc"] in l), None)
                    if i is not None:
                        j = next((j for j in range(i, len(txt)) if txt[j].startswith("Work dir:")), i + 12)
                        fl["err"] = "\n".join(txt[i:j]).strip()
                fails.append(fl)
            if self.status == "failed" and not fails:  # not a task failure (bad path, bad param): show nextflow's own ERROR text
                msgs = list(dict.fromkeys(l[0].strip() for l in self.log if l[0].startswith("ERROR ~")))
                msgs = [m for m in msgs if "pipeline failed" not in m] or msgs
                if msgs:
                    fails = [{"proc": "nextflow", "workdir": "", "err": "\n".join(msgs), "cmd": self.cmd}]
            return {
                "elapsed": round(time.time() - self.t0, 1) if self.t0 else 0,
                "outdir": self.outdir, "status": self.status, "code": self.code, "cmd": self.cmd,
                "rows": rows, "failures": fails, "done": sum(r["done"] for r in rows),
                "total": sum(r["total"] for r in rows),
                "log": self.log[since:], "next": len(self.log),
            }


RUN = Run()


class Handler(BaseHTTPRequestHandler):
    def _send(self, body, ctype="application/json", code=200):
        data = body.encode() if isinstance(body, str) else body
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _local_host(self):
        # Host check blocks DNS-rebinding reads of the directory listing.
        return re.match(r"^(127\.0\.0\.1|localhost)(:\d+)?$", self.headers.get("Host", ""))

    def do_GET(self):
        path, _, query = self.path.partition("?")
        if not self._local_host():
            return self._send("forbidden", "text/plain", 403)
        if path == "/":
            page = (HTML.read_text().replace("__FIELDS__", json.dumps(FIELDS))
                    .replace("__BROWSE__", json.dumps(BROWSE))
                    .replace("__PROFILES__", json.dumps(PROFILES))
                    .replace("__ENV__", json.dumps({
                        "nxf_ver": os.environ.get("NXF_VER", NXF_VER),
                        "java_home": os.environ.get("JAVA_HOME", JAVA_HOME)})))
            self._send(page, "text/html; charset=utf-8")
        elif path == "/api/state":
            m = re.search(r"since=(\d+)", query)
            self._send(json.dumps(RUN.snapshot(int(m.group(1)) if m else 0)))
        elif path == "/api/results":
            self._send(json.dumps(results(parse_qs(query).get("outdir", [""])[0])))
        elif path == "/api/ls":
            self._send(json.dumps(list_dir(parse_qs(query).get("path", [""])[0])))
        else:
            self._send("not found", "text/plain", 404)

    def do_POST(self):
        # Reject cross-site POSTs: any web page could otherwise hit localhost
        # and launch a pipeline run.
        origin = self.headers.get("Origin")
        if not self._local_host() or (
            origin and not re.match(r"^https?://(127\.0\.0\.1|localhost)(:\d+)?$", origin)
        ):
            return self._send('{"errors":["forbidden"]}', code=403)
        n = int(self.headers.get("Content-Length", 0))
        try:
            params = json.loads(self.rfile.read(n) or b"{}")
        except ValueError:
            return self._send('{"errors":["bad JSON"]}', code=400)
        if self.path == "/api/preview":
            args, errors = build_command(params)
            self._send(json.dumps({"cmd": shlex.join(args) if args else "", "errors": errors}))
        elif self.path == "/api/run":
            errors = RUN.start(params)
            self._send(json.dumps({"errors": errors}), code=400 if errors else 200)
        elif self.path == "/api/open":
            ok = open_file(params.get("outdir", ""), params.get("file", ""))
            self._send(json.dumps({"ok": ok}))
        elif self.path == "/api/stop":
            RUN.stop()
            self._send("{}")
        else:
            self._send("not found", "text/plain", 404)

    def log_message(self, *a):
        pass


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--no-browser", action="store_true")
    ap.add_argument("--replay", metavar="LOG", help="show a saved nextflow log as a finished run (parser/UI testing)")
    a = ap.parse_args()
    if a.replay:
        for ln in open(a.replay, errors="replace"):
            RUN._line(ln)
        RUN.status, RUN.code, RUN.t0 = "done", 0, time.time() - 3600
    srv = ThreadingHTTPServer(("127.0.0.1", a.port), Handler)
    url = f"http://localhost:{a.port}"
    print(f"bagRNA GUI on {url}  (Ctrl+C to quit)")
    if not a.no_browser:
        threading.Timer(0.5, webbrowser.open, args=(url,)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
