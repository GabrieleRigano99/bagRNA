#!/bin/bash
# Versions required by bagRNA (see ../CLAUDE.md "Running the pipeline").
export NXF_VER=25.04.6
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
exec python3 "$(dirname "$(readlink -f "$0")")/bagrna_web.py" "$@"
