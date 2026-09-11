#!/usr/bin/env bash

# Pull all required Docker images for the bagRNA pipeline

echo "Pulling Docker images..."

docker pull ftricomi/cpc2:latest
docker pull gabrielerigano/effectorp3:latest
docker pull quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0
docker pull interpro/interproscan:5.75-106.0
docker pull quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_2
docker pull antismash/standalone:8.0.4
docker pull nextgenusfs/funannotate:v1.8.17
docker pull ezlabgva/busco:v6.0.0_cv1
docker pull baderlab/mikado:ubuntu22_mikado2.3.2
docker pull quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2
echo "✅ All images pulled successfully."
