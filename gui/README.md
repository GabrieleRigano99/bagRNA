# bagRNA web GUI

Browser front end for bagRNA. One stdlib Python server (`bagrna_web.py`) plus one page
(`bagrna_web.html`). No Flask, no npm, no display needed. It sits beside `../launch_gui.py`.

## Start

```bash
gui/start_bagrna.sh                 # sets NXF_VER=25.04.6 + JAVA_HOME (java 21), opens the browser
gui/start_bagrna.sh --no-browser --port 8765
ssh -L 8765:localhost:8765 host     # remote use: forward the port, then open http://localhost:8765
```

The server binds `127.0.0.1` only and rejects foreign `Host` / `Origin` headers, because it launches processes.

## What it does

- **Mode** selector: *Run annotation* or *Set up databases* (`-entry SETUP`, shows only `--db_dir` + `--skip_*_db`).
- **Form** is generated at startup from `nextflow.config` `params {}` (defaults, section groups) decorated with
  `nextflow_schema.json` (descriptions, enums, file/dir pickers). A new param shows up with no GUI edit.
  The schema lags the pipeline, so `nextflow.config` wins where they differ.
- **Command** box shows the exact `nextflow run main.nf ...` line, rebuilt on every edit. Only values that differ
  from the default are sent. Validation copies `check_required()` in `main.nf`; Nextflow stays the authority.
- **Run / Stop**: Stop sends SIGTERM to Nextflow's whole process group. The run starts in the pipeline directory,
  so relative paths (e.g. `--outdir run_X`) resolve there.
- **Progress**: bar, process table, per-stage flow diagram (click a chip to filter the log), failure panel
  (stderr tail + `.command.sh` of the failed task).
- **Results**: *Load results* reads `<outdir>/annotation_report/*_annotation_summary.tbl` of any finished run;
  *Open report* / *Open output folder* use `xdg-open`.

## Where things come from

| Thing | Source |
|---|---|
| Process names (un-ellipsizing `STR…`) | `process` / `include … as` in `modules/`, `workflows/`, `subworkflows/` + static `tag "…"` |
| Stage of each process | `include` blocks of `workflows/*.nf`; `subworkflows/` = InterProScan 6 |
| Status lines | `STATUS_RE` / `PEND_RE` in `bagrna_web.py` (Nextflow's progress table). Edit these if Nextflow's output format changes |

Known limits: rows are keyed by process name, so the three processes that exist twice inside InterProScan 6
(`RUN_PFSEARCH`, `PARSE_PFSEARCH`, `REPORT_NO_MATCHES`) share one row. One run at a time. No RNA-seq coverage tile
(not in the summary table).

`--replay LOG` shows a saved Nextflow log as a finished run (parser / UI testing).
