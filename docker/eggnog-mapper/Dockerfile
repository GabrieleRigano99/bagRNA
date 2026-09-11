# eggNOG-mapper 3.0.0-beta5 (eggNOG DB v7), built from the same
# apptainer/eggnog-mapper.def used to produce the upstream .sif.
#
# Build context must be the eggnog-mapper git checkout itself (not this
# directory), e.g.:
#   docker build -f docker/eggnog-mapper/Dockerfile \
#       -t gabrielerigano/eggnog-mapper:3.0.0-beta5 \
#       /mnt/newvolume/new_egg7/eggnog-mapper

FROM python:3.12-slim-bookworm

ARG DIAMOND_VERSION=2.0.15
ARG MMSEQS_VERSION=latest
ARG HMMER_VERSION=3.4
ARG PRODIGAL_VERSION=2.6.3

ENV PATH=/opt/diamond:/opt/mmseqs/bin:/opt/hmmer/bin:/opt/prodigal:$PATH \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
        ca-certificates curl gcc g++ make zlib1g-dev libgomp1 \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/diamond \
    && curl -sL https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VERSION}/diamond-linux64.tar.gz | tar -xz -C /opt/diamond \
    && chmod +x /opt/diamond/diamond

RUN mkdir -p /opt/mmseqs \
    && if [ "$MMSEQS_VERSION" = "latest" ]; then \
        curl -sL https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz | tar -xz -C /opt; \
    else \
        curl -sL https://github.com/soedinglab/MMseqs2/releases/download/${MMSEQS_VERSION}/mmseqs-linux-avx2.tar.gz | tar -xz -C /opt; \
    fi

RUN mkdir -p /tmp/hmmer-build \
    && curl -sL http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz | tar -xz -C /tmp/hmmer-build --strip-components=1 \
    && cd /tmp/hmmer-build \
    && ./configure --prefix=/opt/hmmer > /dev/null \
    && make -j$(nproc) > /dev/null \
    && make install > /dev/null \
    && cd / && rm -rf /tmp/hmmer-build

RUN mkdir -p /opt/prodigal \
    && curl -sL https://github.com/hyattpd/Prodigal/releases/download/v${PRODIGAL_VERSION}/prodigal.linux -o /opt/prodigal/prodigal \
    && chmod +x /opt/prodigal/prodigal

RUN python3.12 -m pip install --break-system-packages --no-cache-dir "cython>=3.0"

COPY . /opt/emapper-src
RUN python3.12 -m pip install --break-system-packages --no-cache-dir /opt/emapper-src \
    && rm -rf /opt/emapper-src

RUN emapper.py --version \
    && diamond --version \
    && mmseqs version \
    && prodigal -v 2>&1 | head -2
