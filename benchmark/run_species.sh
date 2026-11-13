#!/usr/bin/env bash
# Run bagRNA (full mode) on one benchmark species, from scratch (no -resume).
# Usage: benchmark/run_species.sh <species> [extra nextflow args, e.g. -preview]
# Edit the PATHS block once for your machine. Inputs come from benchmark/fetch_inputs.sh.
set -euo pipefail
sp=${1:?usage: run_species.sh <species> [nextflow args]}; shift || true

# ---- PATHS (edit) ----
INPUTS=${INPUTS:-/mnt/newvolume/bagRNA_repro/inputs}      # fetch_inputs.sh output
RUNS=${RUNS:-/mnt/newvolume/bagRNA_repro/runs}            # results + work per species
DBS=${DBS:-/mnt/newvolume/download_test}                  # SETUP output (Rfam, dbCAN, MEROPS, ...)
EGGNOG=${EGGNOG:-/mnt/newvolume/new_egg7}                 # eggNOG v7 data dir
IPS6=${IPS6:-/mnt/newvolume/databases_bag/interproscan-5.78-109.0/data/}
SIGNALP=${SIGNALP:-/home/astrogab/tools/signalp-6.0i.fast/signalp6_fast/}   # licensed, manual download
PHOBIUS=${PHOBIUS:-/home/astrogab/tools/phobius/}                           # licensed, manual download
SBT=${SBT:-/home/astrogab/tools/submission_template.sbt}  # NCBI template made on the submission portal
THREADS=${THREADS:-20}
# ----------------------
# species | strain | locus_tag | BUSCO lineage
case $sp in
  Aspergillus_aculeatinus)   strain=CBS121060; tag=ASPACU; busco=aspergillus ;;
  Candida_orthopsilosis)     strain=90-125;    tag=CANORT; busco=debaryomycetaceae ;;
  Cryptococcus_depauperatus) strain=CBS7841;   tag=CRYDEP; busco=tremellomycetes ;;
  Malassezia_globosa)        strain=CBS7966;   tag=MALGLO; busco=basidiomycota ;;
  Sporothrix_schenckii)      strain=1099-18;   tag=SS109918; busco=sordariomycetes ;;
  Syncephalastrum_racemosum) strain=NRRL2496;  tag=SYNRAC; busco=mucorales ;;
  *) echo "unknown species $sp" >&2; exit 1 ;;
esac
# InterProScan6 merge steps run inside the Nextflow JVM; give it room (default heap = 1/4 of RAM)
export NXF_OPTS=${NXF_OPTS:--Xmx24g}
export NXF_VER=${NXF_VER:-25.04.6} JAVA_HOME=${JAVA_HOME:-/usr/lib/jvm/java-21-openjdk-amd64}
here=$(cd "$(dirname "$0")/.." && pwd); i=$INPUTS/$sp; o=$RUNS/$sp; mkdir -p "$o"
cd "$here"
nextflow run main.nf "$@" \
  --genome_fasta $i/genome.fna --prot_evidence $i/prot_evidence.faa --star_manifest $i/star_manifest.tsv \
  --busco_lineage $busco --species $sp --strain $strain --locus_tag $tag \
  --submission_template $SBT --databases $DBS --dbcan_db $DBS/dbcan_db \
  --merops_db $DBS/merops/merops_pepunit.dmnd --phi_base_db $DBS/phi_base/phi_base.dmnd \
  --go_obo /mnt/newvolume/databases_bag/fun_db/go.obo \
  --gene2product /mnt/newvolume/databases_bag/fun_db/ncbi_cleaned_gene_products.txt \
  --eggnog_data_dir $EGGNOG --IPS6_databases_path $IPS6 --interproscan6_interpro_version 109.0 \
  --signalp_path $SIGNALP --phobius_path $PHOBIUS \
  --helixer_lineage fungi --annevo_lineage fungi --use_gpu --extra_docker_mounts /mnt/newvolume \
  --threads $THREADS --outdir $o/results -w $o/work
