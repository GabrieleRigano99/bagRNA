# bagRNA web GUI

A browser front end for the bagRNA pipeline: fill in parameters, run or stop the pipeline, watch live
per-process progress, see why a task failed, and load results from a finished run.

- One stdlib Python server (`bagrna_web.py`) plus one static page (`bagrna_web.html`).
- No Flask, no npm, no display needed (works over `ssh -L`).
- The GUI never re-implements pipeline logic. It builds a `nextflow run main.nf ...` command, runs it,
  and reads Nextflow's own terminal output.

## Quick start

```bash
cd gui
./start_bagrna.sh                  # sets NXF_VER=25.04.6 and JAVA_HOME (Java 21), opens the browser
./start_bagrna.sh --no-browser     # same, without opening a browser window
./start_bagrna.sh --port 9000      # default port is 8765
```

Machine-specific settings: put `export NXF_WORK=/path/to/work` (the Nextflow work directory; `-resume` only
finds a cache in the same one) in `gui/local.env`. It is sourced by `start_bagrna.sh`, ignored by git, and shows up
in the "NXF_WORK" field of the launch-environment box. Blank means `./work`.

Remote machine: `ssh -L 8765:localhost:8765 <host>`, start the GUI there, open http://localhost:8765 locally.

Requirements: Python 3 (stdlib only), Nextflow, Java 21, and Docker. Tick `-stub-run` in the GUI to try it
without running any real analysis.

## How it works

```
 Browser (bagrna_web.html)                    bagrna_web.py (127.0.0.1 only)            Pipeline
 ┌──────────────────────────┐   POST /api/run  ┌──────────────────────────┐  spawn  ┌───────────────┐
 │ 1. Mode + parameter form │ ───────────────► │ build_command()          │ ──────► │ nextflow run  │
 │ 2. Command preview       │ ◄─────────────── │ validate (as main.nf)    │         │ main.nf ...   │
 │ 3. Run / Stop            │   /api/preview   └──────────────────────────┘         └───────┬───────┘
 │                          │                                                               │ stdout
 │ 4. Progress bar, table,  │   GET /api/state ┌──────────────────────────┐                 │
 │    flow diagram, log,    │ ◄─────────────── │ Run: parse Nextflow's    │ ◄───────────────┘
 │    failure panel         │   every ~1 s     │ progress table + errors  │
 │                          │                  └──────────────────────────┘
 │ 5. Results dashboard     │   GET /api/results ┌────────────────────────┐
 │                          │ ◄───────────────── │ read finished output   │ ◄── <outdir>/annotation_report/
 └──────────────────────────┘                    └────────────────────────┘
```

The page polls `GET /api/state?since=N` about once a second and receives only the new log lines plus the process
table. Stop sends SIGTERM to Nextflow's whole process group, so its containers stop too.

## Using it

| Step | What to do |
|---|---|
| 1. Mode | **Run annotation**, **Structural annotation only** (a normal run with `--no_functional_anno`; hides the functional-database inputs), or **Set up databases** (`-entry SETUP`, shows only `--db_dir` and the `--skip_*_db` flags). |
| 2. Fill the form | Required fields first. **Browse** opens a server-side file picker. Advanced groups are behind a toggle. |
| 3. Check the command | The box shows the exact command, rebuilt on every edit. Errors for missing required fields appear above it. |
| 4. Run | Tick `-stub-run` for a dry end-to-end test, `-resume` to reuse cached results. Then press **Run**. |
| 5. Watch | Progress bar and ETA, process table, flow diagram (click a chip to filter the log), raw log. |
| 6. Failures | The panel shows the failed task's `.command.err` tail and `.command.sh`. If there is no such file, it shows Nextflow's own error text. |
| 7. Results | **Load results** reads any finished output directory (not only the current run). **Open report** and **Open output folder** use `xdg-open`. |

Run status is one of `idle`, `running`, `done`, `failed`, `stopped` (you pressed Stop).
The run starts in the pipeline directory, so a relative `--outdir run_X` resolves there.

## Where the data comes from

Nothing is hand-typed, so the GUI follows the pipeline when it changes.

| GUI part | Source |
|---|---|
| Form fields, defaults, groups | `params {}` in `nextflow.config` (wins on any conflict), described by `nextflow_schema.json` |
| Required-field rules | `check_required()` in `main.nf` |
| Process names (also un-ellipsizing `STR…`) | `process` and `include … as` in `modules/`, `workflows/`, `subworkflows/`, plus static `tag "…"` values |
| Stage of each process | `include` blocks of `workflows/*.nf`; `subworkflows/` counts as InterProScan 6 |
| Progress lines | `STATUS_RE` and `PEND_RE` in `bagrna_web.py` (edit these if Nextflow's output format changes) |
| Results tiles | `<outdir>/annotation_report/*_annotation_summary.tbl`; without it, compleasm `summary.txt` and final GFF feature counts |

## Safety

- Binds `127.0.0.1` only. Every request must carry a local `Host` header (DNS-rebinding guard) and every POST a
  local `Origin` (cross-site launch guard).
- Arguments are passed to `nextflow` as a list, never through a shell.
- One run at a time. No authentication, since only the local machine can connect.

## Known limits

- Rows are keyed by process name, so processes that exist twice inside InterProScan 6
  (`RUN_PFSEARCH`, `PARSE_PFSEARCH`, `REPORT_NO_MATCHES`) share one row.
- `-stub-run` cannot pass the InterProScan 6 `SPLIT_FASTA` step, so a stub run of full functional annotation
  fails there. Use `--no_functional_anno` for stub tests.
- No RNA-seq coverage tile (that figure is not in the summary table).
- Slurm and SGE profiles are passed through as `-profile`, but have not been tested on a cluster.

## Testing the parser

`./start_bagrna.sh --replay path/to/run.log` shows a saved Nextflow log as a finished run, so the parser and
page can be tested without running anything.
