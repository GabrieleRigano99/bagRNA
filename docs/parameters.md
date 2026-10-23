# Parameters reference

Generated from `nextflow_schema.json`. Also available via `nextflow run main.nf --help`.

## Mandatory inputs

| Parameter | Type | Description |
|---|---|---|
| `--genome_fasta` | string | Genome assembly FASTA file. |
| `--species` | string | Species name used for annotation (e.g. `Fusarium_oxysporum`). |
| `--prot_evidence` | string | Protein evidence FASTA for Miniprot and DIAMOND. |
| `--busco_lineage` | string | BUSCO lineage name (e.g. `sordariomycetes`, `fungi_odb10`). |
| `--star_manifest` | string | TSV: `R1.fastq.gz<TAB>R2.fastq.gz<TAB>sample_id`. Repeat one line per sample. |
| `--mikado_config` | string | Mikado list TSV: one transcript source per line. |
| `--submission_template` | string | NCBI `.sbt` submission template file. |

## Helixer options

| Parameter | Type | Description |
|---|---|---|
| `--helixer_lineage` | string | Helixer lineage for de novo GPU prediction (`fungi`, `land_plant`, `vertebrate`, `invertebrate`). |
| `--helixer_gff` | string | Precomputed Helixer GFF — skips the GPU run. |
| `--no_helixer` | boolean | Skip Helixer entirely. |

## ANNEVO options

| Parameter | Type | Description |
|---|---|---|
| `--annevo_lineage` | string | ANNEVO lineage model. |
| `--ram_annevo` | string | RAM limit for ANNEVO (e.g. `40gb`). |
| `--use_gpu` | boolean | Enable GPU for Helixer, ANNEVO, TMbed, and SignalP6. |

## Optional inputs

| Parameter | Type | Description |
|---|---|---|
| `--scoring` | string | Mikado scoring YAML (built-in name or path to custom file). Default: `assets/fungi.yaml`. |
| `--transcript_evidence` | string | Transcript evidence FASTA; aligned with GMAP and passed to Mikado prepare. |
| `--lifted_annotation` | string | Liftover annotation GFF (e.g. from Liftoff). |
| `--phobius_path` | string | Path to Phobius install (mounted at `/opt/phobius`). |
| `--databases` | string | Path to databases directory (EggNOG, KEGG, RFAM). |
| `--signalp_path` | string | Path to `signalp-6-package/` (mounted at `/tools`). Omit to skip SignalP6. |

## InterProScan6 options

| Parameter | Type | Description |
|---|---|---|
| `--IPS6_databases_path` | string | Path to local InterProScan data dir (e.g. `interproscan-5.75-106.0/data`). |
| `--interproscan6_interpro_version` | string | InterPro data version — pin to match `databases_path`. |
| `--interproscan6_apps` | string | Comma-separated analyses to run (e.g. `pfam,panther,smart`). Empty = all. |
| `--interproscan6_goterms` | boolean | Include GO term cross-references. |
| `--interproscan6_pathways` | boolean | Include pathway cross-references (MetaCyc, Reactome). |
| `--interproscan6_no_matches_api` | boolean | Disable InterPro Matches API lookup — force all analyses local. |
| `--interproscan6_tmbed_signalp6_gpu_batch_size` | integer | TMbed/SignalP6 GPU batch size (GPU batch = this × 10). Lower if GPU OOMs. |
| `--no_interpro` | boolean | Skip InterProScan6 entirely. |

## Functional-annotation-only mode

| Parameter | Type | Description |
|---|---|---|
| `--functional_anno_only` | boolean | Run only functional annotation (skips structural). |
| `--protein_fasta` | string | Protein FASTA — required with `--functional_anno_only`. |
| `--ncrna_fasta` | string | ncRNA transcript FASTA — optional, used by Infernal. |
| `--final_gff` | string | Structural annotation GFF — optional, passed to funannotate annotate. |
| `--gbk` | string | table2asn `.gbf` file — optional, passed to AntiSMASH. |

## RNA-seq options

| Parameter | Type | Description |
|---|---|---|
| `--orientation` | string | Read orientation: `FR`, `RF`, or `unstranded`. |
| `--strandedness` | string | Library strandedness for StringTie. |

## NCBI submission

| Parameter | Type | Description |
|---|---|---|
| `--strain` | string | Isolate/strain name. |
| `--locus_tag` | string | Locus tag prefix for GFF features. |
| `--codon_table` | integer | NCBI genetic code table number (1 = standard, 4 = Mycoplasma/Spiroplasma, etc.). |

## Performance

| Parameter | Type | Description |
|---|---|---|
| `--threads` | integer | Number of CPU threads. |
| `--max_memory` | string | Cap memory for high-resource processes (e.g. `50GB`). |
| `--ram_trinity` | string | Memory limit for Trinity de novo assembly (e.g. `45G`). |
| `--limitBAMsortRAM` | integer | STAR BAM sorting RAM in bytes. Auto-computed as 90% of task memory if unset. |
| `--max_intron_length` | integer | Max intron length, passed to Mikado and STAR. |
| `--max_gene_length` | integer | Max gene length used in AGAT GFF filtering. |
| `--genomeSAindexNbases` | integer | STAR `genomeSAindexNbases` — reduce for small genomes (max 14). |

## Feature flags

| Parameter | Type | Description |
|---|---|---|
| `--jaccard_clip` | boolean | Enable Jaccard clip in Trinity (recommended for compact genomes). |
| `--no_functional_anno` | boolean | Skip all functional annotation steps. |
| `--no_antismash` | boolean | Skip AntiSMASH secondary metabolite prediction. |
| `--no_effectorp3` | boolean | Skip EffectorP-3 effector prediction. |
| `--no_eggnog` | boolean | Skip EggNOG-mapper annotation. |
| `--no_signalp` | boolean | Skip SignalP6 signal peptide prediction. |

## Output

| Parameter | Type | Description |
|---|---|---|
| `--outdir` | string | Directory where results are published. |
