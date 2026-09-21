// structural_annotation.nf
// Structural annotation sub-workflow: stages 1-40
// Covers genome BUSCO through protein BUSCO QC

nextflow.enable.dsl = 2

include { COMPLEASM_GENOME                } from '../modules/compleasm_genome'
include { MINIMAP2_ALIGN                  } from '../modules/minimap2_align'
include { MINIPROT                        } from '../modules/miniprot'
include { SUBSAMPLE_BAM                   } from '../modules/subsample_bam'
include { INFER_STRANDEDNESS              } from '../modules/infer_strandedness'
include { HELIXER                         } from '../modules/helixer'
include { ANNEVO_GFF2GTF                  } from '../modules/annevo_gff2gtf'
include { STAR_INDEX                      } from '../modules/star_index'
include { STAR_ALIGN                      } from '../modules/star_align'
include { SAMTOOLS_INDEX                  } from '../modules/samtools_index'
include { SAMTOOLS_SPLIT                  } from '../modules/samtools_split'
include { PORTCULLIS                      } from '../modules/portcullis'
include { MERGE_JUNCTIONS                 } from '../modules/merge_junctions'
include { ALETSCH                         } from '../modules/aletsch'
include { STRINGTIE_ASSEMBLE              } from '../modules/stringtie_assemble'
include { TRINITY                         } from '../modules/trinity'
include { GMAP_BUILD                      } from '../modules/gmap_build'
include { GMAP_ALIGN_TRINITY              } from '../modules/gmap_align_trinity'
include { GMAP_ALIGN_TRANSCRIPTS          } from '../modules/gmap_align_transcripts'
include { AGAT_FIX_CDS_LIFTOVER          } from '../modules/agat_fix_cds_liftover'
include { CORRECT_MISORIENTED             } from '../modules/correct_misoriented'
include { FILTER_BY_JUNCTIONS             } from '../modules/filter_by_junctions'
include { REMOVE_READTHROUGHS             } from '../modules/remove_readthroughs'
include { FILTER_ISOFORMS                 } from '../modules/filter_isoforms'
include { FILTER_UNSUPPORTED              } from '../modules/filter_unsupported'
include { TD2_GENOME_GFF                  } from '../modules/td2_genome_gff'
include { ADD_ISOFORMS                    } from '../modules/add_isoforms'
include { CLEAN_ISOFORMS_GFF              } from '../modules/clean_isoforms_gff'
include { FINAL_JUNCTION_FILTER           } from '../modules/final_junction_filter'
include { MIKADO_PREPARE                  } from '../modules/mikado_prepare'
include { DIAMOND_MAKEDB                  } from '../modules/diamond_makedb'
include { DIAMOND_BLASTX                  } from '../modules/diamond_blastx'
include { TRANSDECODER2_LONGORFS          } from '../modules/transdecoder2'
include { TRANSDECODER2_PREDICT           } from '../modules/transdecoder2'
include { TRANSDECODER2_PFAM_SCAN         } from '../modules/transdecoder2_pfam'
include { TRANSDECODER2_BLASTP            } from '../modules/transdecoder2_blastp'
include { CPC2                            } from '../modules/cpc2'
include { SALMON_INDEX                    } from '../modules/salmon_index'
include { SALMON_QUANT                    } from '../modules/salmon_quant'
include { COMPLEASM_TRANSCRIPTS           } from '../modules/compleasm_transcripts'
include { BUILD_EXTERNAL_SCORES           } from '../modules/build_external_scores'
include { MIKADO_SERIALISE                } from '../modules/mikado_serialise'
include { MIKADO_PICK                     } from '../modules/mikado_pick'
include { FILTER_CODING_MODELS            } from '../modules/filter_coding_models'
include { EXTRACT_NONCODING_MODELS        } from '../modules/extract_noncoding_models'
include { ANNEVO                          } from '../modules/annevo'
include { BARRNAP                         } from '../modules/barrnap'
include { TRNASCAN                        } from '../modules/trnascan'
include { AGAT_FIX_OVERLAPS              } from '../modules/agat_fix_overlaps'
include { MERGE_ADDITIONAL_MODELS         } from '../modules/merge_additional_models'
include { GFF_AGAT_FILTER                 } from '../modules/gff_agat_filter'
include { ADD_MISSING_LOCI                } from '../modules/add_missing_loci'
include { COMPLEASM_MISSING_BUSCO         } from '../modules/compleasm_missing_busco'
include { ADD_BUSCO_ISOFORMS              } from '../modules/add_busco_isoforms'
include { FIX_MICRO_INTRONS               } from '../modules/fix_micro_introns'
include { FUNANNOTATE_RENAME              } from '../modules/funannotate_rename'
include { AGAT_RENAME_IDS as AGAT_RENAME_IDS_1                   } from '../modules/agat_rename_ids'
include { AGAT_RENAME_IDS as AGAT_RENAME_IDS_2                   } from '../modules/agat_rename_ids'
include { REFORMAT_LOCUS_TAG_IDS as REFORMAT_LOCUS_TAG_IDS_1     } from '../modules/reformat_locus_tag_ids'
include { REFORMAT_LOCUS_TAG_IDS as REFORMAT_LOCUS_TAG_IDS_2     } from '../modules/reformat_locus_tag_ids'
include { GFF_CLEAN_FILTER                } from '../modules/gff_clean_filter'
include { TABLE2ASN                       } from '../modules/table2asn'
include { AGAT_EXTRACT_PROTEINS           } from '../modules/agat_extract_proteins'
include { AGAT_EXTRACT_PROTEINS as AGAT_EXTRACT_PROTEINS_INTERIM } from '../modules/agat_extract_proteins'
include { AGAT_EXTRACT_NCRNA              } from '../modules/agat_extract_ncrna'
include { AGAT_MERGE_ABINITIO              } from '../modules/agat_merge_abinitio'
include { GFFCOMPARE_DEDUP_MINIPROT        } from '../modules/gffcompare_dedup_miniprot'

workflow STRUCTURAL_ANNOTATION {

    take:
    ch_fasta                // tuple val(meta), path(fasta)
    ch_prot_evidence        // path
    ch_samples              // channel of tuples: val(meta), path(r1), path(r2)  — used by Salmon
    ch_star_manifest        // collected path(s): manifest file(s) passed to STAR --readFilesManifest
    ch_scoring              // val: scoring filename (built-in) or basename of custom file
    ch_scoring_file         // path: custom scoring file, or NO_SCORING sentinel
    ch_transcript_evidence  // path: transcript evidence FASTA, or NO_FILE sentinel
    ch_submission_template  // path
    ch_lr_reads             // collected FASTQ paths for long reads, or NO_FILE sentinel
    busco_lineage           // val
    species                 // val
    strain                  // val

    main:

    // ── 1. Genome compleasm + miniprot annotation ───────────────────────────
    COMPLEASM_GENOME(ch_fasta, busco_lineage)

    // ── 2. Helixer gene prediction (optional, used as additional models) ──────
    if (params.no_helixer) {
        ch_helixer_gff = Channel.value([ [id: 'genome'], file("${projectDir}/assets/NO_FILE") ])
    } else if (params.helixer_gff) {
        ch_helixer_gff = ch_fasta.map { meta, fasta ->
            [ meta, file(params.helixer_gff) ]
        }
    } else if (params.helixer_lineage) {
        HELIXER(ch_fasta, params.helixer_lineage, species)
        ch_helixer_gff = HELIXER.out.helixer_gff
    } else {
        log.error "Helixer is required. Provide '--helixer_gff' (precomputed) or '--helixer_lineage' (de novo GPU run), or skip with '--no_helixer'."
        System.exit(1)
    }

    // ── 2b. ANNEVO gene prediction → GTF (primary source for STAR/Mikado) ────
    ch_fasta_plain = ch_fasta.map { meta, fasta -> fasta }

    MINIPROT(ch_fasta_plain, ch_prot_evidence)

    ANNEVO(ch_fasta_plain, species, strain, params.annevo_lineage)

    ANNEVO_GFF2GTF(
        ANNEVO.out.annevo_gff.map { gff -> [ [id: 'genome'], gff ] }
    )
    ch_annevo_gtf = ANNEVO_GFF2GTF.out.annevo_gtf

    // ── 2c. Barrnap rRNA + tRNAscan-SE tRNA annotation ────────────────────────
    // Genome-only (no RNA-seq needed) — run unconditionally here so both the
    // full RNA-seq path and the ab-initio-only (no --star_manifest) path
    // below can use them without running twice.
    BARRNAP(ch_fasta)
    AGAT_FIX_OVERLAPS(BARRNAP.out.barrnap_gff)

    TRNASCAN(ch_fasta)

    ch_gtf_plain = ch_annevo_gtf.map { meta, gtf -> gtf }

    // ── 3+. RNA-seq-driven structural annotation (Mikado consensus), or ──────
    //       ab-initio-only fallback when no --star_manifest is given ─────────
    if (params.star_manifest) {

    // ── 3. STAR genome index ────────────────────────────────────────────────
    STAR_INDEX(ch_fasta_plain, ch_gtf_plain)

    // ── 4. STAR alignment (all samples via manifest, RG tag per sample) ─────
    // Collect all read files so Nextflow mounts their directories inside Docker
    ch_star_reads = ch_samples
        .flatMap { meta, r1, r2 -> r2 ? [r1, r2] : [r1] }
        .collect()

    STAR_ALIGN(ch_star_manifest, STAR_INDEX.out.index, ch_star_reads)

    // ── 4b. Index merged STAR BAM for Portcullis ─────────────────────────────
    ch_bam_indexed = SAMTOOLS_INDEX(
        STAR_ALIGN.out.bam.map { bam -> [ [id: 'star'], bam ] }
    )

    // ── 4c. Subsample BAM and infer strandedness ──────────────────────────────
    SUBSAMPLE_BAM(ch_bam_indexed.bam_bai)

    INFER_STRANDEDNESS(
        SUBSAMPLE_BAM.out.bam_bai,
        ch_annevo_gtf.map { meta, gtf -> gtf }
    )

    // User can override auto-detection with --strandedness / --orientation
    if (params.strandedness) {
        ch_strandedness = Channel.value(params.strandedness)
        ch_orientation  = Channel.value(params.orientation ?: 'FR')
        def _tlib = params.strandedness == 'firststrand' ? 'RF'
                  : params.strandedness == 'secondstrand' ? 'FR' : ''
        ch_trinity_lib = Channel.value(_tlib)
        log.info "Strandedness overridden by user: strandedness=${params.strandedness}  orientation=${params.orientation ?: 'FR'}"
    } else {
        ch_strandedness = INFER_STRANDEDNESS.out.strandedness_file.map    { it.text.trim() }
        ch_orientation  = INFER_STRANDEDNESS.out.orientation_file.map     { it.text.trim() }
        ch_trinity_lib  = INFER_STRANDEDNESS.out.trinity_lib_type_file.map { it.text.trim() }
    }

    // ── 5. Portcullis splice site filtering (on merged BAM) ──────────────────
    PORTCULLIS(
        ch_bam_indexed.bam_bai,
        ch_fasta_plain,
        ch_orientation,
        ch_strandedness
    )

    // ── 6. Merge junctions ──────────────────────────────────────────────────
    ch_junction_beds = PORTCULLIS.out.junctions
        .map { meta, bed -> bed }
        .collect()

    MERGE_JUNCTIONS(ch_junction_beds)

    // ── 7. Split Portcullis filtered BAM by RG (condition) → per-condition BAMs
    SAMTOOLS_SPLIT(PORTCULLIS.out.filtered_bam)

    ch_filtered_bams = SAMTOOLS_SPLIT.out.split_bams
        .flatten()
        .collect()

    // ── 7b. Long-read alignment with minimap2 (optional) ─────────────────────
    if (params.lr_manifest) {
        MINIMAP2_ALIGN(
            ch_fasta_plain,
            ch_lr_reads,
            MERGE_JUNCTIONS.out.merged_bed
        )
        ch_lr_bam = MINIMAP2_ALIGN.out.bam
    } else {
        ch_lr_bam = Channel.value(file("${projectDir}/assets/NO_FILE"))
    }

    // ── 8a. Aletsch transcript assembly (per-condition filtered BAMs) ─────────
    ALETSCH(ch_filtered_bams)

    // ── 8b. StringTie transcript assembly; --mix with long reads when available
    STRINGTIE_ASSEMBLE(
        PORTCULLIS.out.filtered_bam.map { meta, bam -> bam },
        ch_gtf_plain,
        ch_strandedness,
        ch_lr_bam
    )

    // ── 8c. Trinity genome-guided assembly (Portcullis filtered BAM) ─────────
    // ── 9. GMAP index, Trinity alignment, and optional transcript evidence ────
    GMAP_BUILD(ch_fasta)
    ch_gmap_index = GMAP_BUILD.out.gmap_index.map { m, i -> i }

    if (!params.no_trinity) {
        TRINITY(
            PORTCULLIS.out.filtered_bam.map { meta, bam -> bam },
            params.jaccard_clip,
            params.ram_trinity,
            ch_trinity_lib
        )
        GMAP_ALIGN_TRINITY(TRINITY.out.trinity_fasta, ch_gmap_index)
        ch_trinity_gff = GMAP_ALIGN_TRINITY.out.trinity_gff
    } else {
        ch_trinity_gff = Channel.value(file("${projectDir}/assets/NO_FILE"))
    }

    if (params.transcript_evidence) {
        GMAP_ALIGN_TRANSCRIPTS(ch_transcript_evidence, ch_gmap_index)
        ch_transcript_evidence_gff = GMAP_ALIGN_TRANSCRIPTS.out.transcript_evidence_gff
    } else {
        ch_transcript_evidence_gff = Channel.value(file("${projectDir}/assets/NO_FILE"))
    }

    // ── 10. Optional liftover annotation ────────────────────────────────────
    if (params.lifted_annotation) {
        AGAT_FIX_CDS_LIFTOVER(
            file(params.lifted_annotation),
            ch_fasta_plain
        )
        ch_liftover_gff = AGAT_FIX_CDS_LIFTOVER.out.liftover_gff
    } else {
        ch_liftover_gff = Channel.value(file("${projectDir}/assets/NO_FILE"))
    }

    // ── 10b. Correct misoriented monoexonic transcripts ──────────────────────
    CORRECT_MISORIENTED(
        STRINGTIE_ASSEMBLE.out.stringtie_gtf,
        ALETSCH.out.aletsch_gtf,
        ch_trinity_gff,
        MINIPROT.out.miniprot_gtf
    )

    // ── 10c. Remove assembled transcripts with Portcullis-unvalidated junctions
    FILTER_BY_JUNCTIONS(
        CORRECT_MISORIENTED.out.stringtie_gtf,
        CORRECT_MISORIENTED.out.aletsch_gtf,
        CORRECT_MISORIENTED.out.trinity_gff,
        PORTCULLIS.out.junctions.map { meta, bed -> bed }
    )

    // ── 11. Mikado prepare ──────────────────────────────────────────────────
    // Combine prot_evidence with BUSCO refseq database sequences
    MIKADO_PREPARE(
        ch_fasta_plain,
        MERGE_JUNCTIONS.out.merged_bed,
        ch_prot_evidence,
        ch_scoring,
        ch_scoring_file,
        ch_helixer_gff.map { meta, gff -> gff },
        ch_annevo_gtf.map   { meta, gtf -> gtf },
        FILTER_BY_JUNCTIONS.out.stringtie_gtf,
        FILTER_BY_JUNCTIONS.out.trinity_gff,
        FILTER_BY_JUNCTIONS.out.aletsch_gtf,
        COMPLEASM_GENOME.out.busco_anno_gff.map { meta, gff -> gff },
        MINIPROT.out.miniprot_gtf,
        ch_liftover_gff,
        ch_transcript_evidence_gff
    )

    // ── 12. Diamond protein database + blastx ───────────────────────────────
    DIAMOND_MAKEDB(ch_prot_evidence)
    DIAMOND_BLASTX(
        MIKADO_PREPARE.out.prepared_fasta,
        DIAMOND_MAKEDB.out.diamond_db
    )

    // ── 13. TD2 ORF prediction with homology-guided retention ────────────────
    TRANSDECODER2_LONGORFS(MIKADO_PREPARE.out.prepared_fasta, ch_strandedness)

    // Pfam hmmscan (if Pfam HMM is available via --pfam_db, --IPS6_databases_path,
    // or a databases/interproscan/pfam dir downloaded by -entry SETUP)
    def pfam_hmm_file = null
    if (params.pfam_db) {
        pfam_hmm_file = file(params.pfam_db)
    } else {
        def pfamBase = params.IPS6_databases_path ?: (params.databases ? "${params.databases}/interproscan" : null)
        if (pfamBase) {
            def pfamDir = file("${pfamBase}/pfam").listFiles()?.find { it.isDirectory() }
            if (pfamDir) pfam_hmm_file = file("${pfamDir}/pfam_a.hmm")
        }
    }

    if (pfam_hmm_file?.exists()) {
        TRANSDECODER2_PFAM_SCAN(TRANSDECODER2_LONGORFS.out.pep, pfam_hmm_file)
        ch_td2_hmmer = TRANSDECODER2_PFAM_SCAN.out.domtblout
    } else {
        ch_td2_hmmer = Channel.value(file("${projectDir}/assets/NO_FILE"))
    }

    // blastp against protein evidence (prot_evidence is already required for
    // structural annotation; a separate eggnog_proteins.dmnd pass used to be
    // appended here too, but it was always supplementary — dropped 2026-08-31).
    TRANSDECODER2_BLASTP(
        TRANSDECODER2_LONGORFS.out.pep,
        DIAMOND_MAKEDB.out.diamond_db
    )

    TRANSDECODER2_PREDICT(
        MIKADO_PREPARE.out.prepared_fasta,
        TRANSDECODER2_LONGORFS.out.td2_dir,
        ch_td2_hmmer,
        TRANSDECODER2_BLASTP.out.blastp_hits
    )

    // ── 14. CPC2 coding potential ────────────────────────────────────────────
    CPC2(MIKADO_PREPARE.out.prepared_fasta)

    // ── 15. Salmon quantification ────────────────────────────────────────────
    SALMON_INDEX(MIKADO_PREPARE.out.prepared_fasta)

    // Collect all R1 and R2 reads from the manifest into single files
    ch_reads_r1 = ch_samples.map { meta, r1, r2 -> r1 }.collect()
    ch_reads_r2 = ch_samples.map { meta, r1, r2 -> r2 }.collect()

    SALMON_QUANT(
        SALMON_INDEX.out.salmon_index,
        ch_reads_r1,
        ch_reads_r2
    )

    // ── 16. External scores (CPC2 + Salmon TPM + BUSCO protein) ─────────────
    // Scan TD2 proteins against BUSCO HMMs using the DB downloaded by COMPLEASM_GENOME
    COMPLEASM_TRANSCRIPTS(
        TRANSDECODER2_LONGORFS.out.pep,
        busco_lineage,
        COMPLEASM_GENOME.out.busco_db
    )

    BUILD_EXTERNAL_SCORES(
        SALMON_QUANT.out.quant_sf,
        CPC2.out.cpc2_result,
        COMPLEASM_TRANSCRIPTS.out.busco_scores
    )

    // ── 17. Mikado serialise + pick ───────────────────────────────────────────
    MIKADO_SERIALISE(
        MIKADO_PREPARE.out.mikado_dir,
        MIKADO_PREPARE.out.configuration,
        ch_scoring,
        ch_scoring_file,
        ch_fasta_plain,
        MERGE_JUNCTIONS.out.merged_bed,
        TRANSDECODER2_PREDICT.out.bed,
        DIAMOND_BLASTX.out.diamond_tsv,
        BUILD_EXTERNAL_SCORES.out.external_scores,
        ch_prot_evidence
    )

    MIKADO_PICK(MIKADO_SERIALISE.out.mikado_dir_serialised, MIKADO_PREPARE.out.configuration, ch_scoring, ch_scoring_file, ch_fasta_plain)

    // ── 17b. Remove readthrough transcripts ──────────────────────────────────
    REMOVE_READTHROUGHS(MIKADO_PICK.out.pick_gff, MINIPROT.out.miniprot_gtf)

    // ── 17c. Keep genuine AS isoforms (primary=False ccode=j), drop the rest ──
    FILTER_ISOFORMS(REMOVE_READTHROUGHS.out.filtered_gff)

    // ── 17c2. Correct multi-exonic ab initio primaries with no portcullis-confirmed junction ──
    FILTER_UNSUPPORTED(
        FILTER_ISOFORMS.out.filtered_gff,
        ch_junction_beds
    )

    // ── 17c3. Lift TD2 transcript-space ORFs to genome coordinates ───────────
    TD2_GENOME_GFF(
        TRANSDECODER2_PREDICT.out.gff3,
        MIKADO_PREPARE.out.prepared_gtf,
        MIKADO_PREPARE.out.prepared_fasta
    )

    // ── 17d. Add TD2-predicted isoforms (all sources) with RNA-seq-verified junctions ──
    ADD_ISOFORMS(
        FILTER_UNSUPPORTED.out.filtered_gff,
        TD2_GENOME_GFF.out.genome_gff3,
        ch_junction_beds,
        SALMON_QUANT.out.quant_sf
    )

    // ── 17e. Fix gene boundaries + deduplicate isoforms ──────────────────────
    CLEAN_ISOFORMS_GFF(ADD_ISOFORMS.out.augmented_gff)

    // ── 17f. Remove any remaining models with unconfirmed junctions ───────────
    FINAL_JUNCTION_FILTER(
        CLEAN_ISOFORMS_GFF.out.cleaned_gff,
        ch_junction_beds
    )

    // ── 18. Post-Mikado processing ────────────────────────────────────────────
    // Reverted from AGAT_RENAME_IDS back to FUNANNOTATE_RENAME (2026-08-28):
    // the AGAT swap (2026-08-26) is directly correlated with a severe bug where
    // AGAT_RENAME_IDS_2 (much further downstream, in the locus-tag rename pass)
    // merges hundreds of unrelated genes into single multi-Mb records — every
    // pre-swap run (mid/late July) is clean, every post-swap run (26+28 Aug)
    // has it, nothing else in the pipeline changed. Root cause inside
    // agat_sp_manage_IDs.pl itself not isolated; reverting this early pass is
    // the well-evidenced fix. Confirmed funannotate's renamer never touches
    // tRNA at this position anyway (tRNA enters later via TRNASCAN), so this
    // doesn't reintroduce the tRNA product=None problem the swap targeted —
    // that protection lives in AGAT_RENAME_IDS_1/_2, unchanged either way.
    FILTER_CODING_MODELS(FINAL_JUNCTION_FILTER.out.filtered_gff)
    FUNANNOTATE_RENAME(FILTER_CODING_MODELS.out.coding_gff, ch_fasta_plain)
    EXTRACT_NONCODING_MODELS(FINAL_JUNCTION_FILTER.out.filtered_gff)

    // ── 22. Merge additional gene models ──────────────────────────────────────
    ch_liftover_input = params.lifted_annotation ? ch_liftover_gff : Channel.value(file("${projectDir}/assets/NO_FILE"))

    MERGE_ADDITIONAL_MODELS(
        ANNEVO.out.annevo_gff,
        AGAT_FIX_OVERLAPS.out.barrnap_gff.map { meta, gff -> gff },
        TRNASCAN.out.trnascan_gff.map { meta, gff -> gff },
        EXTRACT_NONCODING_MODELS.out.noncoding_gff,
        FUNANNOTATE_RENAME.out.renamed_gff,
        ch_helixer_gff.map { meta, gff -> gff },
        ch_liftover_input
    )

    // ── 23. GFF clean, filter, rename ─────────────────────────────────────────
    GFF_AGAT_FILTER(
        MERGE_ADDITIONAL_MODELS.out.merged_gff,
        ch_fasta_plain,
        params.max_gene_length
    )

    // ── 22b. Recover real BUSCO genes lost/mis-modeled by evidence gaps or
    //         AGAT's incomplete-ORF filtering — same approach validated
    //         manually on SS02 (see HANDOFF.md), now run on every genome.
    //         Recovered loci still pass through GFF_CLEAN_FILTER's NCBI QC
    //         and locus-tag renaming below like everything else.
    ADD_MISSING_LOCI(
        GFF_AGAT_FILTER.out.filtered_gff,
        TD2_GENOME_GFF.out.genome_gff3,
        COMPLEASM_GENOME.out.busco_anno_gff.map { meta, gff -> gff }
    )

    AGAT_EXTRACT_PROTEINS_INTERIM(
        ADD_MISSING_LOCI.out.augmented_gff,
        ch_fasta_plain,
        params.codon_table
    )

    COMPLEASM_MISSING_BUSCO(
        AGAT_EXTRACT_PROTEINS_INTERIM.out.proteins_faa,
        busco_lineage,
        COMPLEASM_GENOME.out.busco_db
    )

    ADD_BUSCO_ISOFORMS(
        ADD_MISSING_LOCI.out.augmented_gff,
        COMPLEASM_GENOME.out.busco_anno_gff.map { meta, gff -> gff },
        COMPLEASM_MISSING_BUSCO.out.missing_ids
    )

    // ── 22c. Merge spurious sub-10bp ab-initio "introns" that would
    //         otherwise cost the whole gene at NCBI QC below (frame-safe,
    //         stop-codon-checked — never touches a gap it can't prove safe).
    FIX_MICRO_INTRONS(
        ADD_BUSCO_ISOFORMS.out.recovered_gff,
        ch_fasta_plain
    )

    ch_pre_rename_gff = FIX_MICRO_INTRONS.out.fixed_gff

    } else {

    // ── Ab-initio-only structural annotation (no --star_manifest) ────────────
    // Skips the entire RNA-seq-driven evidence chain (STAR, Portcullis,
    // Aletsch/StringTie/Trinity, Mikado, TransDecoder2, Salmon) — there is no
    // transcript evidence to build any of it from. Gene models come straight
    // from Helixer + ANNEVO (ab initio) and Miniprot (protein-to-genome);
    // rRNA/tRNA from Barrnap/tRNAscan-SE; merged and cleaned with AGAT
    // instead of Mikado's evidence-weighted consensus. No isoform recovery,
    // no BUSCO-gap recovery (both depend on TransDecoder2/RNA-seq junctions).
    log.warn "No '--star_manifest' provided: running ab-initio-only structural annotation " +
             "(Helixer + ANNEVO + Miniprot + Barrnap + tRNAscan-SE, merged with AGAT). " +
             "No transcript assembly, no Mikado consensus, no isoform/BUSCO-gap recovery."

    GFFCOMPARE_DEDUP_MINIPROT(
        ch_helixer_gff.map { meta, gff -> gff },
        ch_annevo_gtf.map  { meta, gtf -> gtf },
        MINIPROT.out.miniprot_gtf,
        AGAT_FIX_OVERLAPS.out.barrnap_gff.map { meta, gff -> gff }
    )

    AGAT_MERGE_ABINITIO(
        ch_helixer_gff.map { meta, gff -> gff },
        ch_annevo_gtf.map  { meta, gtf -> gtf },
        GFFCOMPARE_DEDUP_MINIPROT.out.miniprot_dedup_gtf,
        AGAT_FIX_OVERLAPS.out.barrnap_gff.map { meta, gff -> gff },
        TRNASCAN.out.trnascan_gff.map { meta, gff -> gff }
    )

    GFF_AGAT_FILTER(
        AGAT_MERGE_ABINITIO.out.merged_gff,
        ch_fasta_plain,
        params.max_gene_length
    )

    ch_pre_rename_gff = GFF_AGAT_FILTER.out.filtered_gff

    }

    // ── 23. Locus-tag rename + NCBI QC (AGAT preserves product=/isotype=/
    //         anticodon=/Name= attributes that funannotate gff-rename used to
    //         destroy — see modules/reformat_locus_tag_ids.py). Renamed twice:
    //         once before tbl2gbk QC, once after gene removal to close the
    //         numbering gaps it leaves — same two-pass shape funannotate used.
    AGAT_RENAME_IDS_1(ch_pre_rename_gff, params.locus_tag)
    REFORMAT_LOCUS_TAG_IDS_1(
        AGAT_RENAME_IDS_1.out.agat_renamed_gff,
        params.locus_tag,
        'renamed_pass1.gff3'
    )

    GFF_CLEAN_FILTER(
        REFORMAT_LOCUS_TAG_IDS_1.out.renamed_gff,
        ch_fasta_plain,
        species,
        strain,
        params.locus_tag
    )

    AGAT_RENAME_IDS_2(GFF_CLEAN_FILTER.out.filtered_gff, params.locus_tag)
    REFORMAT_LOCUS_TAG_IDS_2(
        AGAT_RENAME_IDS_2.out.agat_renamed_gff,
        params.locus_tag,
        "final_struct_${species}_${strain}.gff"
    )

    // ── 24. Table2asn (NCBI submission files) ─────────────────────────────────
    TABLE2ASN(
        REFORMAT_LOCUS_TAG_IDS_2.out.renamed_gff,
        ch_fasta_plain,
        ch_submission_template,
        species,
        strain,
        params.codon_table,
        params.locus_tag
    )

    // ── 25. Extract final proteins and ncRNA transcripts ──────────────────────
    AGAT_EXTRACT_PROTEINS(
        REFORMAT_LOCUS_TAG_IDS_2.out.renamed_gff,
        ch_fasta_plain,
        params.codon_table
    )

    AGAT_EXTRACT_NCRNA(
        REFORMAT_LOCUS_TAG_IDS_2.out.renamed_gff,
        ch_fasta_plain
    )

    emit:
    final_gff          = REFORMAT_LOCUS_TAG_IDS_2.out.renamed_gff
    proteins_faa       = AGAT_EXTRACT_PROTEINS.out.proteins_faa
    ncrna_fasta        = AGAT_EXTRACT_NCRNA.out.ncrna_fasta
    gbk                = TABLE2ASN.out.gbf_files
    busco_genome_summary = COMPLEASM_GENOME.out.summary
}
