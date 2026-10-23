#!/bin/bash
# Versions required by bagRNA (see ../CLAUDE.md "Running the pipeline").
export NXF_VER=25.04.6
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
# Machine-specific settings (e.g. export NXF_WORK=/big/disk/work) go in local.env, which git ignores.
here="$(dirname "$(readlink -f "$0")")"
[ -f "$here/local.env" ] && . "$here/local.env"
exec python3 "$here/bagrna_web.py" "$@"
