# Structural annotation

`workflows/structural_annotation.nf` — genome + RNA-seq evidence → gene models.

## 1. Genome QC

`COMPLEASM_GENOME` — BUSCO genome-mode completeness check.

## 2. Ab initio predictions

- `HELIXER` — GPU, optional (`--no_helixer` to skip, `--helixer_gff` to supply a precomputed one)
- `ANNEVO` — GPU; its GTF is used to build the primary STAR index

## 3. Protein-to-genome alignment

`MINIPROT`, using `--prot_evidence`.

## 4. RNA-seq alignment

`STAR_INDEX` → `STAR_ALIGN`. All samples in `--star_manifest` are aligned together in one STAR run via `--readFilesManifest`, with per-sample RG tags.

## 5. Splice-junction filtering

`PORTCULLIS` → `MERGE_JUNCTIONS`. A downsampled BAM is used only for strandedness inference; the full BAM goes to Portcullis.

## 6. Transcript assembly

- `ALETSCH` — per-condition assembly
- `STRINGTIE_ASSEMBLE` — with optional long-read mixing (`--mix`)
- `TRINITY` — genome-guided

## 7. Evidence integration

`MIKADO_PREPARE` collects all GFFs/GTFs into one list, then:

- `DIAMOND_BLASTX`
- `TRANSDECODER2_LONGORFS` / `TRANSDECODER2_PREDICT` — GPU, uses PSAURON neural network for ORF scoring
- `CPC2`
- `SALMON_QUANT`
- `BUILD_EXTERNAL_SCORES`

## 8. Gene model selection

`MIKADO_SERIALISE` → `MIKADO_PICK`, scored against `assets/fungi.yaml` (or `--scoring`). Fungal-tuned vs stock Mikado defaults: lowered `combined_cds_fraction` (`>0.2`), `as_requirements` enabled (keeps isoforms), 500 bp UTR filters, rescaled `exon_num`.

## 9. Post-Mikado processing

- `FILTER_CODING_MODELS` — drops `ncRNA` and `mRNA` models missing a start/stop codon (`has_start_codon=False` / `has_stop_codon=False`). This is a filter, not a repair — flagged models are removed, not fixed.
- `FUNANNOTATE_RENAME` — renames gene/transcript IDs to locus tags
- `BARRNAP` — rRNA detection
- `TRNASCAN` — tRNA detection
- `MERGE_ADDITIONAL_MODELS` — merges Mikado models with BARRNAP/tRNAscan models

## 10. Final cleanup

`GFF_AGAT_FILTER` → `GFF_CLEAN_FILTER` → `TABLE2ASN` → `AGAT_EXTRACT_PROTEINS` / `AGAT_EXTRACT_NCRNA`.

## Other GFF-touching steps

| Process | What it fixes |
|---|---|
| `AGAT_FIX_CDS_LIFTOVER` | CDS coordinates after a liftover |
| `FILTER_BY_JUNCTIONS` | drops transcripts unsupported by splice junctions |
| `FILTER_ISOFORMS` | filters isoforms |
| `FILTER_UNSUPPORTED` | removes unsupported models |
| `CLEAN_ISOFORMS_GFF` | cleans isoform GFF |
| `FINAL_JUNCTION_FILTER` | last-pass junction filtering |
| `AGAT_FIX_OVERLAPS` | overlapping BARRNAP rRNA features |
| `FIX_MICRO_INTRONS` | spuriously short introns |
| `AGAT_RENAME_IDS` (×2) | ID renaming passes |
| `AGAT_MERGE_ABINITIO` | merges ab initio predictor models into final GFF |

## Known issue

**BUSCO 622109 dropout** — gene `mikado.CBS145945_7G858` (7 exons, ANNEVO-based) is correctly picked by Mikado but silently removed between Mikado pick and the final GFF. Under investigation.
