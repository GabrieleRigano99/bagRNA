# Functional annotation

`workflows/functional_annotation.nf` — gene models → functional evidence + submission-ready output.

Runs automatically after structural annotation in full mode, or standalone with `--functional_anno_only`.

## Processes

| Process | What it does | Skip flag |
|---|---|---|
| `PREPARE_INTERPROSCAN` | stages input for the IPS6 subworkflow (`cache = false` — input is not Kryo-serializable) | `--no_interpro` |
| IPS6 subworkflow | InterProScan6 across member databases: Pfam, CDD, PIRSR, MobiDBLite, SignalP, TMBed, DeepTMHMM, etc. — one container per member DB | `--no_interpro` |
| `EGGNOG` | eggNOG-mapper functional orthology | `--no_eggnog` |
| `DEEPKOALA` | KEGG Orthology (KO) assignment — GRU-based, replaced KofamScan (60–107× faster, higher recall) | — |
| `KEGG_ANNOTATE` / `KEGG_KO_INFO` | reformats DeepKOALA output to KofamScan's TSV layout, resolves KO metadata | — |
| `FETCH_KEGG_MODULE_COMPLETENESS` | KEGG module completeness for the report | — |
| `INFERNAL` | Rfam ncRNA family search | — |
| `SIGNALP6` | signal peptide prediction (needs `--signalp_path`) | `--no_signalp` |
| `PHOBIUS` | transmembrane/signal peptide prediction (needs `--phobius_path`) | — |
| `RUN_DBCAN` | CAZyme annotation | — |
| `MEROPS` | peptidase family annotation | — |
| `PHI_BASE` | pathogen-host interaction database matching | — |
| `EFFECTORP3` | fungal effector prediction | `--no_effectorp3` |
| `ANTISMASH` | secondary metabolite biosynthetic gene cluster prediction (needs `--gbk`) | `--no_antismash` |
| `ANNOTATE_FUNCTIONAL` | merges all functional evidence into the final annotation | — |
| `GFF3_TO_TBL` | converts to NCBI `.tbl` format | — |
| `GENERATE_REPORT` | self-contained HTML + PDF annotation report (`gabrielerigano/bagrna-report`, chromium baked in) | — |

`--no_functional_anno` skips this whole workflow.

## DeepKOALA notes

Replaced KofamScan 2026-09-14 (benchmarked 60–107× faster, higher recall at comparable precision, cross-validated against eggNOG's independent KO calls on 2 species). Runs fine on CPU — GPU access is opportunistic (`--gpus all` only if `--use_gpu` is set elsewhere), unlike the pipeline's other `use_gpu`-labeled steps. Its optional `--multi` domain-validation mode (not currently wired in) would need the KOfam HMM database re-added.

## AntiSMASH known issue

`--functional_anno_only` without `--gbk` warns it will skip AntiSMASH but currently does not. Workaround: pass `--no_antismash` explicitly.

