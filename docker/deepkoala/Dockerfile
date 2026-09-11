# DeepKOALA — deep learning-based KEGG Orthology (KO) assignment
# (https://github.com/zhaoxi120/deepkoala)
#
# Build context must be the deepkoala git checkout itself (not this
# directory), e.g.:
#   docker build -f docker/deepkoala/Dockerfile \
#       -t gabrielerigano/deepkoala:0.1-beta \
#       /path/to/deepkoala

FROM pytorch/pytorch:2.4.1-cuda12.1-cudnn9-runtime

ARG HMMER_VERSION=3.4

ENV PATH=/opt/hmmer/bin:$PATH \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# hmmer (hmmsearch) is only needed for the optional --multi domain-validation
# mode; the KOfam profile database itself is not bundled here — mount it at
# runtime and pass its path via --profiles_dir.
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
        ca-certificates curl gcc make zlib1g-dev libgomp1 \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /tmp/hmmer-build \
    && curl -sL http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz | tar -xz -C /tmp/hmmer-build --strip-components=1 \
    && cd /tmp/hmmer-build \
    && ./configure --prefix=/opt/hmmer > /dev/null \
    && make -j$(nproc) > /dev/null \
    && make install > /dev/null \
    && cd / && rm -rf /tmp/hmmer-build

# torch is already provided by the base image (2.4.1+cu121, satisfies
# requirements.txt's torch>=2.4.1,<3); strip it out before installing the
# rest so pip doesn't reinstall/upgrade it.
COPY requirements.txt /tmp/requirements.txt
RUN grep -v '^torch' /tmp/requirements.txt > /tmp/requirements.notorch.txt \
    && pip install --no-cache-dir -r /tmp/requirements.notorch.txt \
    && rm /tmp/requirements.txt /tmp/requirements.notorch.txt

# Editable install: keeps deepkoala/ (the package) and resources/ (the
# pretrained weights) as siblings on disk, exactly matching infer.py's
# Path(__file__).resolve().parent.parent / "resources" lookup. A normal
# `pip install .` would copy only the package into site-packages and
# separate it from resources/, breaking model-weight resolution.
WORKDIR /opt/deepkoala
COPY . /opt/deepkoala
RUN pip install --no-cache-dir --no-deps -e .

RUN deepkoala --help > /dev/null \
    && hmmsearch -h > /dev/null

# No ENTRYPOINT: Nextflow always invokes containers as
# `docker run <image> /bin/bash -c "<generated .command.run>"` — a fixed
# ENTRYPOINT would swallow that whole invocation as arguments to itself
# instead of running as a shell (confirmed the hard way: DEEPKOALA's first
# real in-pipeline run failed with deepkoala's own argparse complaining
# -i/-o were "missing" — they were present but never reached because
# `/bin/bash -c "..."` itself became deepkoala's argv). Matches
# docker/eggnog-mapper and docker/td2, neither of which sets one either;
# call the `deepkoala` console script explicitly instead (see module).
