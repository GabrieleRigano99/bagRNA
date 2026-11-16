# bagRNA v1.0.0 reproducibility report

Date: 2026-11-14. Run by Gabriele Rigano (BCGL-Unime) with Claude Code as assistant.

The question: can someone who starts from an empty folder, only this repository and the public
sources, rebuild the results for all six benchmark fungi, and do the runs finish cleanly?
Each species was run from scratch (no `-resume`) in full mode. The pass rule is below.

## Result

All six species pass. Each run ended with exit code 0 and no failed task, and the NCBI
`table2asn` validation file has no `Error:` line for any species.

| Species | Commit | Wall time | Tasks | NCBI errors / warnings | Genes (Sep) | BUSCO protein % (Sep) | % genes with KO (Sep) | % genes with InterPro (Sep) |
|---|---|---|---|---|---|---|---|---|
| Malassezia globosa | `378e6d9` | 4 h 37 min | 361 | 0 / 1,268 | 4,942 (4,904) | 82.61 (82.40) | 70.4 (70.2) | 80.5 (80.2) |
| Candida orthopsilosis | `ba628b7` (1) | 5 h 54 min | 346 | 0 / 299 | 6,365 (6,358) | 99.04 (99.04) | 67.5 (67.5) | 79.3 (79.3) |
| Cryptococcus depauperatus | `378e6d9` | 18 h 33 min | 758 | 0 / 4,443 | 8,417 (8,192) | 96.60 (94.42) | 52.3 (52.4) | 63.0 (63.0) |
| Aspergillus aculeatinus | `378e6d9` | 10 h 13 min | 758 | 0 / 2,314 | 17,381 (17,434) | 98.91 (98.59) | 39.6 (39.3) | 57.6 (57.2) |
| Sporothrix schenckii | `378e6d9` | 12 h 52 min | 657 | 0 / 1,497 | 17,209 (17,123) | 96.30 (96.08) | 32.4 (31.7) | 44.9 (44.0) |
| Syncephalastrum racemosum | `34310bf` (2) | 11 h 18 min | 934 | 0 / 3,538 | 14,022 (13,900) | 98.32 (97.19) | 51.2 (51.0) | 65.5 (65.3) |

1. Candida was run on `ba628b7`. `378e6d9` differs from it only in the time-limit lines of
   `nextflow.config` (limits now scale with the attempt number). I checked this with a diff and did not rerun
   Candida, which would cost about 6 h.
2. Syncephalastrum failed on `378e6d9` (see flag 7 below) and was rerun from scratch on `34310bf`.
   `34310bf` differs from `378e6d9` only in a concurrency limit and the Nextflow heap size. The other five species
   were not rerun on it.

"Sep" is the first bagRNA run of each species in September 2026. Totals: 63 h of wall time over six runs.

### Pass rule

For each species:

1. Exit code 0, every task `COMPLETED` in `pipeline_info/execution_trace.txt`.
2. Zero `Error:` lines in the `table2asn` `.val` file. Warnings are allowed (see below).
3. Gene count within 10 % of the September baseline.
4. % of genes with a KEGG KO, and % with an InterPro hit, each within 2 points of the baseline.
5. Protein BUSCO (`compleasm protein`) not more than 1 point below the baseline.

Rule 5 is one-sided. I first wrote it as "within 1 point" in both directions. Cryptococcus
(+2.18) and Syncephalastrum (+1.13) are higher than September, so a two-sided rule would fail them
for doing better.
`benchmark/check_outputs.py` applies the rule as written here.

### About the warnings

Nearly all warnings are intron and exon shape notices: non-consensus splice donor and acceptor,
short exon, and the same checks on partial features. They are expected for gene models predicted from RNA-seq in fungi
and do not block submission. Cryptococcus also has 14 `BadEcNumberValue` warnings.
I accepted three EC-value warnings; I did not fix them.

## Problems found and fixed

Each was found by running the pipeline, then fixed and rerun. The full history, with the commit
for each, is in `REPRO_LOG.md`.

1. **EC number on "hypothetical protein"** (Malassezia, 11 errors). `table2asn -c fx` renames products after the
   table is written. Fix: `bin/strip_ec_from_renamed.sh` and a second `table2asn` pass.
2. **`GeneXrefStrandProblem`** (Candida 4, Cryptococcus 13, Aspergillus 3). One `/gene` symbol sat on many loci.
   Fix: `bin/gff3_to_tbl.py` writes `/gene` only for symbols that belong to one locus. Shared symbols go in the
   CDS note.
3. **`FeatureBeginsOrEndsInGap`** (Candida 2, Aspergillus 2). Fix: `gff3_to_tbl.py` leaves out genes with a boundary
   inside a run of 10 or more N.
4. **Internal stop codons** in BUSCO-recovery (`compleasm`) models, which funannotate `tbl2gbk` does not flag
   (Candida 1 model, Cryptococcus 3). Fix: `bin/drop_bad_isoforms.py --all-internal-stops`, called from
   `modules/gff_clean_filter.nf`.
5. **`NoStop`** (Sporothrix). Fix: a CDS whose length is not a multiple of 3 is now marked 3' partial.
6. **TMBED hit its fixed 2 h limit** under heavy machine load and failed Malassezia after 739 min.
   Fix: the time labels in `nextflow.config` scale with the attempt number.
7. **Java out-of-memory in Nextflow** (Syncephalastrum, after 12 h). Three InterProScan6 steps
   (`REPRESENTATIVE_LOCATIONS_LOCAL`, `XREFS_LOCAL`, `COMBINE_MATCHES_LOCAL`) run inside the Nextflow JVM and each
   holds a file of about 1 GB (JSON) in memory, about 6 GB. Several at once overflowed the default heap (about 15 GB on a 62 GB machine).
   Fix: `maxForks = 1` for these steps in `nextflow.config`, and `NXF_OPTS=-Xmx24g` by default in
   `benchmark/run_species.sh`. A single task alone ran in 17 s with a 6 GB heap. The rerun passed.
8. **Documentation and scripts:** the README dry run now uses `-preview` (`-stub-run` fails at the
   InterProScan6 `SPLIT_FASTA` step); `fetch_inputs.sh` uses NCBI SRA instead of ENA (ENA gave about 110 KB/s here) and
   keeps no partial UniProt file after a failed request.

## What you need

- Linux, Docker with the NVIDIA Container Toolkit, Java 17 or newer, Nextflow 25.04.x, `git`, `curl`, `unzip`,
  `sra-tools`, conda with `compleasm` for the BUSCO check.
- The machine used here: Intel i9-13900HX (32 threads), 62 GB RAM, one RTX 4090 Laptop GPU (16 GB), Docker 27.5.0,
  Nextflow 25.04.6, Java 21. Disk: one `work/` folder was 47 GB (Aspergillus) to 62 GB (Syncephalastrum), and the downloaded SRA files
  were about 45 GB for all six species. I did not measure the total peak.
- Reference databases from the `-entry SETUP` workflow, InterProScan 5.78-109.0 data, eggNOG v7, and
  the two licensed tools SignalP 6 and Phobius (manual download, see the README).
- An NCBI submission template (`.sbt`) made on the submission portal.

## Exact commands

```bash
git clone https://github.com/GabrieleRigano99/bagRNA.git && cd bagRNA
git checkout <commit from the table>          # see README: git "dubious ownership" note

benchmark/fetch_inputs.sh /data/inputs                 # genomes, UniProt evidence, SRA reads, SHA256SUMS
INPUTS=/data/inputs RUNS=/data/runs benchmark/run_species.sh Malassezia_globosa   # one species, from scratch

# BUSCO on the final proteins, then the pass rule
compleasm protein -p /data/runs/Malassezia_globosa/results/structural_annotation/final_proteins.faa \
  -l basidiomycota -L <busco_downloads> -t 6 -o Malassezia_globosa     # run in /data/busco
python3 benchmark/check_outputs.py /data/runs Malassezia_globosa --busco /data/busco
```

The paths at the top of `benchmark/run_species.sh` (databases, licensed tools, template) must be set once for the
machine. `benchmark/run_all.sh` queues several species one after another.
BUSCO lineages: Aspergillus `aspergillus`, Candida `debaryomycetaceae`, Cryptococcus `tremellomycetes`,
Malassezia `basidiomycota`, Sporothrix `sordariomycetes`, Syncephalastrum `mucorales`.

## Known limits

1. **The database download was not tested.** The `-entry SETUP` workflow was skipped on purpose. The EBI mirror gave
   50 to 100 KB/s here. The runs used databases that already existed on the machine, so the download path a new
   user takes is untested.
2. **SignalP 6 and Phobius** need a licence and a manual download.
3. **UniProt evidence** is release 2026_04, not the February 2026 files used in September. Evidence therefore differs
   slightly from the baseline runs.
4. **Sporothrix reads:** September used locally trimmed reads; this run used the raw SRA reads, and results still match.
5. **Candida and Syncephalastrum** were run on commits that differ from `378e6d9` (see table notes), not on one single
   commit. One commit for all six would need roughly 6 h more for Candida and about 11 h more for Syncephalastrum
   plus a rerun of the other four species for the memory fix.
6. **Gene counts are not fully deterministic.** Malassezia gave 4,941 and 4,942 genes in two runs. I did not
   investigate.
7. The `interpro/interpro-n` image does not exist on Docker Hub. InterProScan6 skips it and the runs are not affected.
8. One run per species. I have no estimate of run-to-run variation, except the one-gene difference above.
