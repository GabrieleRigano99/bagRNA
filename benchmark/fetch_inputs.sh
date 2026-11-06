#!/usr/bin/env bash
# Fetch all inputs for the 6-species bagRNA benchmark. Needs: bash, curl, unzip, gzip, sra-tools (prefetch, fasterq-dump).
# Usage: benchmark/fetch_inputs.sh /path/to/inputs [species_regex]   (default: all 6 species)
# Writes <dir>/<species>/{genome.fna, prot_evidence.faa, reads/*.fastq.gz, star_manifest.tsv}
set -euo pipefail
OUT=${1:?usage: fetch_inputs.sh OUT_DIR [species_regex]}; ONLY=${2:-.}
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)

# species | genome accession | UniProt taxon id | RNA-seq runs (name:run, ...)
# Evidence taxa: ascomycota 4890, basidiomycota 5204, mucoromycota 1913637, fungi 4751.
# Aspergillus uses the ascomycota set, Cryptococcus uses all fungi (as in the paper runs).
TABLE='Aspergillus_aculeatinus|GCA_003184765.1|4890|Acu:SRR1799574
Candida_orthopsilosis|GCA_000315875.1|4890|Cortho1:SRR10251160,Cortho2:SRR10251161,Cortho3:SRR10251162
Cryptococcus_depauperatus|GCA_001720195.2|4751|Cdep1:SRR17871193,Cdep2:SRR17871194
Malassezia_globosa|GCA_000181695.2|5204|Mglob:SRR6152146
Sporothrix_schenckii|GCF_000961545.1|4890|yeast:SRR8986163,yeast:SRR8986160,yeast:SRR8986159,mold:SRR8986162,mold:SRR8986161,mold:SRR8986164
Syncephalastrum_racemosum|GCA_002105135.1|1913637|Srac:SRR8843209'

get() { curl -fL --retry 5 --retry-delay 10 --retry-all-errors -sS "$@"; }

echo "$TABLE" | grep -E "^($ONLY)" | while IFS='|' read -r sp acc tax runs; do
  d="$OUT/$sp"; mkdir -p "$d/reads"; echo "== $sp"
  if [ ! -s "$d/genome.fna" ]; then
    get -o "$d/g.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA"
    unzip -p "$d/g.zip" "ncbi_dataset/data/$acc/*_genomic.fna" > "$d/genome.fna"; rm "$d/g.zip"
  fi
  if [ ! -s "$d/prot_evidence.faa" ]; then
    # UniProt sometimes answers 429/502 under load: retry the whole request, never keep a partial file
    for try in 1 2 3 4 5 6; do
      get -o "$d/prot_evidence.tmp" "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28taxonomy_id%3A$tax%29+AND+%28reviewed%3Atrue%29" && [ -s "$d/prot_evidence.tmp" ] && break
      echo "UniProt request failed (try $try), retrying in 60 s" >&2; sleep 60
    done
    [ -s "$d/prot_evidence.tmp" ] || { echo "ERROR: UniProt download failed for taxon $tax" >&2; exit 1; }
    mv "$d/prot_evidence.tmp" "$d/prot_evidence.faa"
    # UniProt release changes over time: record it
    get -sI "https://rest.uniprot.org/uniprotkb/search?query=accession:P69905&size=1" | grep -i '^x-uniprot-release:' > "$d/uniprot_release.txt" || true
  fi
  : > "$d/star_manifest.tsv"; n=0
  for item in ${runs//,/ }; do
    cond=${item%%:*}; run=${item##*:}; n=$((n+1))
    # NCBI SRA (S3 mirror, ~25x faster than ENA FASTQ here). Needs sra-tools (prefetch, fasterq-dump).
    r1="$d/reads/${run}_1.fastq.gz"; r2="$d/reads/${run}_2.fastq.gz"
    if [ ! -s "$r2" ]; then
      prefetch --max-size 100G -O "$d/sra" "$run" >&2
      fasterq-dump --split-files --threads 8 -O "$d/reads" "$d/sra/$run" >&2
      [ -s "$d/reads/${run}_2.fastq" ] || { echo "ERROR: $run is not paired" >&2; exit 1; }
      gzip -f "$d/reads/${run}_1.fastq" "$d/reads/${run}_2.fastq"
    fi
    printf '%s\t%s\t%s\n' "$r1" "$r2" "$cond" >> "$d/star_manifest.tsv"
  done
  (cd "$d" && sha256sum genome.fna prot_evidence.faa reads/*.fastq.gz > SHA256SUMS)
done
echo "DONE. Inputs in $OUT"
