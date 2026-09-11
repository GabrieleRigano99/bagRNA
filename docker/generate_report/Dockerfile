# bagRNA annotation report generator — python:3.11-slim + Chromium pre-baked
# in at build time (as root, during `docker build`).
#
# bin/generate_annotation_report.py used to `apt-get install chromium` at
# container RUNTIME, which worked in every manual test but never in the real
# pipeline: Nextflow always runs containers as the host's non-root UID
# (`-u $(id -u):$(id -g)`), so apt-get always failed with "Permission
# denied" on /var/lib/apt/lists — silently, since the PDF render is
# best-effort (see render_pdf() in the script). Confirmed on both the
# Cryptococcus depauperatus and Syncephalastrum racemosum runs, 2026-09-09:
# HTML/.tbl wrote fine, PDF silently never appeared. Baking chromium in here
# avoids the runtime-root problem entirely.
#
# Build: docker build -t gabrielerigano/bagrna-report:latest docker/generate_report

FROM python:3.11-slim

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
        chromium ca-certificates \
    && rm -rf /var/lib/apt/lists/*
