# Usage

Every run needs these two variables exported first:

```bash
export NXF_VER=25.04.6
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
```

bagRNA has four run modes.

## 1. Database setup (one-time)

Downloads EggNOG, KEGG, RFAM, and other reference databases used by functional annotation.

```bash
nextflow run main.nf -entry SETUP --db_dir /path/to/databases
```

## 2. Full mode — structural + functional annotation

Genome FASTA in, annotated GFF + protein set + functional report out.

```bash
nextflow run main.nf \
    --genome_fasta maskedSS02.fa \
    --prot_evidence proteins.faa \
    --busco_lineage sordariomycetes \
    --star_manifest STAR_config.tsv \
    --species "Sporothrix_schenckii" \
    --submission_template template.sbt \
    --helixer_gff helixer.gff \
    --databases /path/to/databases \
    --IPS6_databases_path /path/to/ips6 \
    --signalp_path /path/to/signalp-6-package \
    --outdir run_NAME \
    --use_gpu \
    -resume
```

`--star_manifest` is a TSV, one line per sample: `R1.fastq.gz<TAB>R2.fastq.gz<TAB>condition`.

Runs `structural_annotation.nf` (genome → gene models) then `functional_annotation.nf` (gene models → functional annotation).

## 3. Functional annotation only

Skip structural prediction; annotate a genome + protein FASTA you already have.

```bash
nextflow run main.nf \
    --functional_anno_only \
    --genome_fasta genome.fa \
    --protein_fasta final_proteins.faa \
    --species "Genus_species" \
    --submission_template template.sbt \
    --databases /path/to/databases \
    --IPS6_databases_path /path/to/ips6
```

!!! warning "AntiSMASH"
    Without `--gbk` (a table2asn `.gbf` file), this mode warns it will skip AntiSMASH but currently does not skip it correctly. Pass `--no_antismash` explicitly as a workaround until the underlying `.nf` fix lands.

## 4. Resume a run

Add `-resume` to any of the above with the same parameters:

```bash
nextflow run main.nf [same params as before] -resume
```

Nextflow reuses cached results for unchanged process definitions and inputs.

!!! danger "Never edit `.nf` files between resumed runs"
    Changing anything under `workflows/`, `modules/`, or `subworkflows/` invalidates Nextflow's cache fingerprint for every downstream process — a resumed run will silently recompute far more than expected. Put per-process overrides (resources, extra Docker mounts) in `nextflow.config` under a `withName:` block instead.

## Execution profiles

```bash
nextflow run main.nf [params] -profile slurm
nextflow run main.nf [params] -profile sge
```

Default is `-profile standard` (local executor).

## GPU usage

`--use_gpu` enables GPU for Helixer, ANNEVO, and TransDecoder2 (hard-required for those). DeepKOALA (KEGG Orthology) uses GPU opportunistically if available — it runs fine on CPU otherwise.

## Full parameter list

See [Parameters reference](parameters.md), or:

```bash
nextflow run main.nf --help
```
