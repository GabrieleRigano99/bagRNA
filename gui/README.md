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
| 1. Mode | **Run annotation**, or **Set up databases** (`-entry SETUP`, shows only `--db_dir` and the `--skip_*_db` flags). |
| 2. Fill the form | Required fields first. **Browse** opens a server-side file picker. Advanced groups are behind a toggle. |
| 3. Check the command | The box shows the exact command, rebuilt on every edit. Errors for missing required fields appear above it. |
| 4. Run | Tick `-stub-run` for a dry end-to-end test, `-resume` to reuse cached results. Then press **Run**. |
| 5. Watch | Progress bar and ETA, process table, flow diagram (click a chip to filter the log), raw log. |
| 6. Failures | The panel shows the failed task's `.command.err` tail and `.command.sh`. If there is no such file, it shows Nextflow's own error text. |
| 7. Results | **Load results** reads any finished output directory (not only the current run). **Open report** and **Open output folder** use `xdg-open`. |

Run status is one of `idle`, `running`, `done`, `failed`, `stopped` (you pressed Stop).
The run starts in the pipeline directory, so a relative `--outdir run_X` resolves there.

