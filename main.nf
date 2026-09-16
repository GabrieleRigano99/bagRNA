#!/usr/bin/env nextflow
// ============================================================================
//  bagRNA — end-to-end noncoding-aware eukaryotic genome annotation pipeline
//  Nextflow DSL2 conversion
//  Author: Gabriele Rigano
//  Bioinformatics and Computational Genomics LAB (BCGL), University of Messina
// ============================================================================

nextflow.enable.dsl = 2

include { STRUCTURAL_ANNOTATION  } from './workflows/structural_annotation'
include { FUNCTIONAL_ANNOTATION  } from './workflows/functional_annotation'
include { DOWNLOAD_DATABASES     } from './workflows/download_databases'

// ── Help message ─────────────────────────────────────────────────────────────
def helpMessage() {
    def R  = "[0m"       // reset
    def B  = "[1m"       // bold
    def DM = "[2m"       // dim
    def GR = "[1;32m"    // bold green
    def CY = "[1;36m"    // bold cyan
    def YL = "[1;33m"    // bold yellow
    def MG = "[1;35m"    // bold magenta
    def LN = "${CY}${'─' * 72}${R}"   // section divider

    log.info """
${GR}   _                 ____   _   _    _${R}
${GR}  | |__   __ _  __ _|  _ \\ | \\ | |  / \\${R}
${GR}  | '_ \\ / _' |/ _' | |_) ||  \\| | / _ \\${R}
${GR}  | |_) | (_| | (_| |  _ < | |\\  |/ ___ \\${R}
${GR}  |_.__/ \\__,_|\\__, |_| \\_\\|_| \\_/_/   \\_\\${R}
${GR}               |___/${R}

  ${B}bagRNA${R} — end-to-end noncoding-aware eukaryotic genome annotation pipeline
  ${DM}Version 1.0.0  •  Bioinformatics & Computational Genomics LAB, UniME${R}
  ${DM}Main Developer  •  Gabriele Rigano - ORCID https://orcid.org/0009-0008-1928-6789 ${R}

  ${B}USAGE:${R}  nextflow run main.nf [options]

${LN}
  ${CY}${B}FULL MODE${R}  —  structural + functional annotation
${LN}

  ${YL}${B}REQUIRED${R}
    ${B}--genome_fasta${R}          Genome assembly FASTA file
    ${B}--prot_evidence${R}         Protein evidence FASTA file
    ${B}--busco_lineage${R}         BUSCO lineage name  ${DM}(e.g. sordariomycetes)${R}
    ${B}--star_manifest${R}         TSV: R1.fastq.gz\\tR2.fastq.gz\\tsample_id
    ${B}--species${R}               Species name  ${DM}(e.g. "Genus_species")${R}
    ${B}--submission_template${R}   NCBI .sbt submission template
    ${B}--databases${R}             Databases directory  ${DM}(KEGG, RFAM)${R}
    ${B}--eggnog_data_dir${R}       Path to eggNOG-mapper data dir  ${DM}(eggnog.db, eggnog_proteins.dmnd, ...)${R}
    ${B}--IPS6_databases_path${R}   Local InterProScan 6 databases directory
    ${B}--signalp_path${R}          Path to signalp-6-package/  ${DM}(mounted at /tools)${R}
    ${B}--phobius_path${R}          Path to Phobius directory  ${DM}(mounted at /opt/phobius)${R}

  ${YL}${B}HELIXER${R}  ${DM}(one required unless --no_helixer)${R}
    ${B}--helixer_gff${R}           Precomputed Helixer GFF  ${DM}(skips GPU run)${R}
    ${B}--helixer_lineage${R}       Lineage for de novo GPU prediction
                          ${DM}fungi | land_plant | vertebrate | invertebrate${R}
    ${B}--no_helixer${R}            Skip Helixer entirely

  ${YL}${B}ANNEVO${R}
    ${B}--annevo_lineage${R}        Model lineage  ${DM}[${params.annevo_lineage}]${R}
                          ${DM}fungi | land_plant | vertebrate | invertebrate${R}
                          ${DM}magnoliopsida | mammalia | insecta${R}

${LN}
  ${CY}${B}FUNCTIONAL ANNOTATION ONLY MODE${R}
${LN}

    ${B}--functional_anno_only${R}  Run only functional annotation  ${DM}(skips structural)${R}
    ${B}--genome_fasta${R}          Genome assembly FASTA file
    ${B}--protein_fasta${R}         Protein FASTA  ${DM}(replaces structural output)${R}
    ${B}--species${R}               Species name
    ${B}--submission_template${R}   NCBI .sbt submission template
    ${B}--databases${R}             Databases directory  ${DM}(KEGG, RFAM)${R}
    ${B}--eggnog_data_dir${R}       Path to eggNOG-mapper data dir  ${DM}(eggnog.db, eggnog_proteins.dmnd, ...)${R}
    ${B}--IPS6_databases_path${R}   Local InterProScan 6 databases directory
    ${B}--signalp_path${R}          Path to signalp-6-package/  ${DM}(mounted at /tools)${R}
    ${B}--phobius_path${R}          Path to Phobius directory  ${DM}(mounted at /opt/phobius)${R}

  ${YL}${B}OPTIONAL${R}  ${DM}(functional_anno_only)${R}
    ${B}--ncrna_fasta${R}           ncRNA FASTA  ${DM}(for Infernal)${R}
    ${B}--final_gff${R}             Structural annotation GFF  ${DM}(for the functional annotation merge)${R}
    ${B}--gbk${R}                   Table2asn .gbf files  ${DM}(for AntiSMASH; omit to skip)${R}
    ${B}--id_map${R}                old_id\\tnew_id TSV  ${DM}(translates pre-rename tool-output IDs)${R}

${LN}
  ${CY}${B}SHARED OPTIONAL INPUTS${R}
${LN}

    ${B}--scoring${R}               Mikado scoring YAML  ${DM}[scerevisiae.yaml]${R}
    ${B}--transcript_evidence${R}   Transcript FASTA  ${DM}(aligned with GMAP → Mikado prepare)${R}
    ${B}--lifted_annotation${R}     Liftover annotation GFF  ${DM}(e.g. from Liftoff)${R}
    ${B}--lr_manifest${R}           TSV of long-read FASTQ paths  ${DM}(one path per line; enables minimap2 + StringTie --mix)${R}
    ${B}--lr_type${R}               Long-read platform preset  ${DM}[${params.lr_type}]  ont | pacbio_hifi${R}
    ${B}--dbcan_db${R}              dbCAN database directory  ${DM}(CAZy.dmnd, dbCAN.hmm, TF.hmm, STP.hmm …)${R}
    ${B}--merops_db${R}             MEROPS pepunit FASTA or pre-built .dmnd  ${DM}(peptidase annotation)${R}
    ${B}--phi_base_db${R}           PHI-base FASTA or pre-built .dmnd  ${DM}(pathogen-host interaction)${R}

${LN}
  ${CY}${B}PERFORMANCE${R}
${LN}

    ${B}--threads${R}               CPU threads  ${DM}[${params.threads}]${R}
    ${B}--max_memory${R}            Memory cap for high-resource processes  ${DM}(e.g. 50GB)${R}
    ${B}--ram_trinity${R}           Trinity memory limit  ${DM}[${params.ram_trinity}]${R}
    ${B}--ram_annevo${R}            ANNEVO memory limit  ${DM}[${params.ram_annevo}]${R}
    ${B}--max_intron_length${R}     Maximum intron length for Mikado and Portcullis  ${DM}[${params.max_intron_length}]${R}
    ${B}--max_gene_length${R}       Maximum gene length  ${DM}[${params.max_gene_length}]${R}
    ${B}--genomeSAindexNbases${R}   STAR SA index bases  ${DM}[${params.genomeSAindexNbases}]  (reduce for small genomes)${R}
    ${B}--limitBAMsortRAM${R}       STAR BAM sort RAM bytes  ${DM}(auto: 90% of task memory)${R}
    ${B}--orientation${R}           Override inferred read orientation: FR | RF | unstranded  ${DM}[auto-detected]${R}
    ${B}--strandedness${R}          Override inferred strandedness: firststrand | secondstrand | unstranded  ${DM}[auto-detected]${R}

${LN}
  ${CY}${B}NCBI SUBMISSION${R}
${LN}

    ${B}--strain${R}                Isolate / strain name  ${DM}[${params.strain}]${R}
    ${B}--locus_tag${R}             GFF locus tag prefix  ${DM}[${params.locus_tag}]${R}
    ${B}--codon_table${R}           Genetic code table  ${DM}[${params.codon_table}]${R}
    ${B}--outdir${R}                Output directory  ${DM}[${params.outdir}]${R}

${LN}
  ${CY}${B}FEATURE FLAGS${R}
${LN}

    ${B}--use_gpu${R}               Enable GPU for Helixer, ANNEVO and TMbed
    ${B}--jaccard_clip${R}          Enable Jaccard clip in Trinity  ${DM}[${params.jaccard_clip}]${R}
    ${B}--no_helixer${R}            Skip Helixer  ${DM}[${params.no_helixer}]${R}
    ${B}--no_functional_anno${R}    Skip all functional annotation
    ${B}--no_interpro${R}           Skip InterProScan 6
    ${B}--no_eggnog${R}             Skip EggNOG
    ${B}--no_antismash${R}          Skip AntiSMASH
    ${B}--no_effectorp3${R}         Skip EffectorP-3
    ${B}--no_signalp${R}            Skip SignalP6
    ${B}--no_dbcan${R}              Skip run_dbcan CAZyme/CGC annotation
    ${B}--no_merops${R}             Skip MEROPS peptidase search
    ${B}--no_phi_base${R}           Skip PHI-base pathogen-host interaction search

${LN}
  ${CY}${B}INTERPROSCAN 6${R}
${LN}

    ${B}--interproscan6_apps${R}              Comma-separated analyses to run  ${DM}(default: all)${R}
                                      ${DM}e.g. pfam,cdd,smart,signalp,phobius${R}
    ${B}--IPS6_databases_path${R}                  Local IPS6 databases directory
    ${B}--interproscan6_interpro_version${R}  InterPro data version  ${DM}[${params.interproscan6_interpro_version}]${R}
                                      ${DM}Pin to match --IPS6_databases_path data (e.g. 106.0)${R}
    ${B}--interproscan6_goterms${R}           Annotate GO terms  ${DM}[${params.interproscan6_goterms}]${R}
    ${B}--interproscan6_pathways${R}          Annotate pathway cross-references  ${DM}[${params.interproscan6_pathways}]${R}
    ${B}--interproscan6_no_matches_api${R}    Skip API lookup for no-match sequences  ${DM}[${params.interproscan6_no_matches_api}]${R}
    ${B}--interproscan6_tmbed_signalp6_gpu_batch_size${R}  TMBed / SignalP6 GPU batch size  ${DM}[${params.interproscan6_tmbed_signalp6_gpu_batch_size}]${R}
                                           ${DM}Actual batch = value × 10; reduce if GPU OOM${R}

${LN}
  ${MG}${B}EXAMPLES${R}
${LN}

  ${DM}# Full run — precomputed Helixer GFF, IPS6 in local mode, GPU enabled${R}
  nextflow run main.nf \\
      --genome_fasta genome.fa                      \\
      --prot_evidence proteins.fa                   \\
      --busco_lineage sordariomycetes               \\
      --star_manifest manifest.tsv                  \\
      --species "Fusarium_oxysporum"                \\
      --submission_template template.sbt            \\
      --helixer_gff helixer.gff                     \\
      --databases /path/to/databases                \\
      --IPS6_databases_path /path/to/ips6_data      \\
      --signalp_path /path/to/signalp-6-package     \\
      --use_gpu

  ${DM}# Full run — de novo Helixer, IPS6 in API mode, no AntiSMASH${R}
  nextflow run main.nf \\
      --genome_fasta genome.fa                      \\
      --prot_evidence proteins.fa                   \\
      --busco_lineage sordariomycetes               \\
      --star_manifest manifest.tsv                  \\
      --species "Fusarium_oxysporum"                \\
      --submission_template template.sbt            \\
      --helixer_lineage fungi                       \\
      --annevo_lineage fungi                        \\
      --databases /path/to/databases                \\
      --no_antismash

  ${DM}# Functional annotation only${R}
  nextflow run main.nf \\
      --functional_anno_only                        \\
      --genome_fasta genome.fa                      \\
      --protein_fasta proteins.faa                  \\
      --species "Fusarium_oxysporum"                \\
      --submission_template template.sbt            \\
      --databases /path/to/databases                \\
      --IPS6_databases_path /path/to/ips6_data      \\
      --signalp_path /path/to/signalp-6-package

  ${DM}# Database setup (run once before the pipeline)${R}
  nextflow run main.nf -entry SETUP --db_dir /path/to/databases
  ${DM}  # skip individual databases${R}
  nextflow run main.nf -entry SETUP --db_dir /path/to/databases --skip_eggnog_db --skip_rfam_db
  ${DM}  # then run the pipeline pointing to the downloaded databases${R}
  nextflow run main.nf ... \\
      --databases      /path/to/databases           \\
      --dbcan_db       /path/to/databases/dbcan     \\
      --merops_db      /path/to/databases/merops/merops_pepunit.dmnd \\
      --phi_base_db    /path/to/databases/phi_base/phi_base.dmnd

${LN}
  ${CY}${B}SETUP FLAGS${R}  ${DM}(only with -entry SETUP)${R}
${LN}

    ${B}--db_dir${R}                Target directory for all databases  ${DM}(required)${R}
    ${B}--skip_eggnog_db${R}        Skip eggNOG-mapper databases  ${DM}(~15 GB)${R}
    ${B}--skip_rfam_db${R}          Skip Rfam covariance models  ${DM}(~1 GB)${R}
    ${B}--skip_dbcan_db${R}         Skip dbCAN databases  ${DM}(~2 GB)${R}
    ${B}--skip_merops_db${R}        Skip MEROPS pepunit + DIAMOND index  ${DM}(~1 GB)${R}
    ${B}--skip_phi_base_db${R}      Skip PHI-base FASTA + DIAMOND index  ${DM}(~6 MB)${R}
    ${B}--skip_go_obo_db${R}        Skip Gene Ontology go-basic.obo  ${DM}(~7 MB)${R}
    ${B}--skip_interproscan_db${R}  Skip InterProScan 6 member databases  ${DM}(~50 GB)${R}

${LN}
    """.stripIndent()
}

// ── Parameter validation ─────────────────────────────────────────────────────
def check_required(param_name, param_val) {
    if (!param_val) {
        log.error "Parameter '--${param_name}' is required but was not provided."
        log.error "Run with '--help' for usage information."
        System.exit(1)
    }
}

// eggNOG-mapper v7 data dir isn't nested under --databases the way KEGG/RFAM
// are read directly by name — but -entry SETUP does publish it to
// db_dir/eggnog_data, so treat that as a fallback for the purposes of
// deciding whether to warn (the actual fallback is resolved again, the same
// way, inside functional_annotation.nf and passed into EGGNOG explicitly).
def eggnogDataDirAvailable() {
    if (params.eggnog_data_dir) return true
    return params.databases && file("${params.databases}/eggnog_data").exists()
}

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {

    if (params.help) {
        helpMessage()
        System.exit(0)
    }

    // ── Mode-aware validation ─────────────────────────────────────────────────
    check_required('genome_fasta', params.genome_fasta)
    check_required('species',      params.species)

    if (params.functional_anno_only) {
        check_required('protein_fasta',        params.protein_fasta)
        check_required('submission_template',  params.submission_template)
    } else {
        check_required('prot_evidence',       params.prot_evidence)
        check_required('busco_lineage',       params.busco_lineage)
        check_required('star_manifest',       params.star_manifest)
        check_required('submission_template', params.submission_template)

        if (!params.no_helixer && !params.helixer_gff && !params.helixer_lineage) {
            log.error "Helixer input required: provide '--helixer_gff' (precomputed GFF), '--helixer_lineage' (de novo GPU run), or '--no_helixer' to skip."
            log.error "Run with '--help' for usage information."
            System.exit(1)
        }
    }

    // ── Print pipeline header ─────────────────────────────────────────────────
    def R  = "\033[0m"
    def B  = "\033[1m"
    def DM = "\033[2m"
    def GR = "\033[1;32m"

    log.info """
${GR}   _                 ____   _   _    _${R}
${GR}  | |__   __ _  __ _|  _ \\ | \\ | |  / \\${R}
${GR}  | '_ \\ / _' |/ _' | |_) ||  \\| | / _ \\${R}
${GR}  | |_) | (_| | (_| |  _ < | |\\  |/ ___ \\${R}
${GR}  |_.__/ \\__,_|\\__, |_| \\_\\|_| \\_/_/   \\_\\${R}
${GR}               |___/${R}

  ${B}bagRNA${R} — end-to-end noncoding-aware eukaryotic genome annotation pipeline
  ${DM}Version 1.0.0  •  Bioinformatics & Computational Genomics LAB, UniME${R}
  ${DM}Main Developer  •  Gabriele Rigano - ORCID https://orcid.org/0009-0008-1928-6789 ${R}


        mode                : ${params.functional_anno_only ? 'functional_anno_only' : 'full'}
        genome_fasta        : ${params.genome_fasta}
        species             : ${params.species}
        strain              : ${params.strain}
        threads             : ${params.threads}
        outdir              : ${params.outdir}
        """.stripIndent()

    // ── Build shared channels ─────────────────────────────────────────────────
    ch_fasta = Channel
        .fromPath(params.genome_fasta, checkIfExists: true)
        .map { fasta -> [ [id: 'genome', species: params.species], fasta ] }

    ch_submission_template = Channel
        .fromPath(params.submission_template, checkIfExists: true)

    if (params.functional_anno_only) {

        // ── Functional-annotation-only mode ───────────────────────────────────
        ch_proteins_faa = Channel
            .fromPath(params.protein_fasta, checkIfExists: true)

        ch_ncrna_fasta = params.ncrna_fasta
            ? Channel.fromPath(params.ncrna_fasta, checkIfExists: true)
            : Channel.value(file("${projectDir}/assets/NO_FILE"))

        ch_final_gff = params.final_gff
            ? Channel.fromPath(params.final_gff, checkIfExists: true)
            : Channel.value(file("${projectDir}/assets/NO_FILE"))

        ch_gbk = params.gbk
            ? Channel.fromPath(params.gbk, checkIfExists: true)
            : Channel.value(file("${projectDir}/assets/NO_FILE"))

        if (!params.databases) {
            log.warn "No '--databases' path provided. KEGG and RFAM searches will be skipped. (InterProScan 6 still runs via API.)"
        }
        if (!params.no_eggnog && !eggnogDataDirAvailable()) {
            log.warn "No '--eggnog_data_dir' provided. EggNOG will be skipped."
        }
        if (!params.gbk && !params.no_antismash) {
            log.warn "No '--gbk' provided in functional_anno_only mode. AntiSMASH will be skipped."
        }

        FUNCTIONAL_ANNOTATION(
            ch_proteins_faa,
            ch_ncrna_fasta,
            ch_final_gff,
            ch_fasta.map { meta, fasta -> fasta },
            ch_gbk,
            ch_submission_template,
            Channel.value(file("${projectDir}/assets/NO_FILE")), // no BUSCO genome-mode run in functional_anno_only mode
            params.species,
            params.strain,
            params.final_gff as boolean,
            params.ncrna_fasta as boolean
        )

    } else {

        // ── Full mode: structural + optional functional annotation ─────────────
        ch_prot_evidence = Channel
            .fromPath(params.prot_evidence, checkIfExists: true)

        def scoring_path = file(params.scoring)
        if (scoring_path.exists()) {
            ch_scoring      = Channel.value(scoring_path.name)
            ch_scoring_file = Channel.fromPath(params.scoring, checkIfExists: true)
        } else {
            ch_scoring      = Channel.value(params.scoring)
            ch_scoring_file = Channel.value(file("${projectDir}/assets/NO_FILE"))
        }

        def manifest_paths = params.star_manifest instanceof List
            ? params.star_manifest
            : [ params.star_manifest ]

        ch_samples = Channel
            .fromPath(manifest_paths, checkIfExists: true)
            .splitCsv(header: false, sep: '\t', strip: true)
            .map { row ->
                def meta = [ id: row[2] ]
                def r1   = file(row[0], checkIfExists: true)
                def r2   = row.size() > 1 && row[1] ? file(row[1], checkIfExists: true) : []
                return [ meta, r1, r2 ]
            }

        ch_star_manifest = Channel
            .fromPath(manifest_paths, checkIfExists: true)
            .collect()

        ch_transcript_evidence = params.transcript_evidence
            ? Channel.fromPath(params.transcript_evidence, checkIfExists: true)
            : Channel.value(file("${projectDir}/assets/NO_FILE"))

        ch_lr_reads = params.lr_manifest
            ? Channel
                .fromPath(params.lr_manifest, checkIfExists: true)
                .splitCsv(header: false, sep: '\t', strip: true)
                .map { row -> file(row[0], checkIfExists: true) }
                .collect()
            : Channel.value(file("${projectDir}/assets/NO_FILE"))

        STRUCTURAL_ANNOTATION(
            ch_fasta,
            ch_prot_evidence,
            ch_samples,
            ch_star_manifest,
            ch_scoring,
            ch_scoring_file,
            ch_transcript_evidence,
            ch_submission_template,
            ch_lr_reads,
            params.busco_lineage,
            params.species,
            params.strain
        )

        if (!params.no_functional_anno) {
            if (!params.databases) {
                log.warn "Functional annotation requested but '--databases' path not provided. " +
                         "Set '--no_functional_anno' to suppress this warning or provide the databases path."
            }
            if (!params.no_eggnog && !eggnogDataDirAvailable()) {
                log.warn "Functional annotation requested but '--eggnog_data_dir' not provided. " +
                         "EggNOG will be skipped. Set '--no_eggnog' to suppress this warning."
            }

            FUNCTIONAL_ANNOTATION(
                STRUCTURAL_ANNOTATION.out.proteins_faa,
                STRUCTURAL_ANNOTATION.out.ncrna_fasta,
                STRUCTURAL_ANNOTATION.out.final_gff,
                ch_fasta.map { meta, fasta -> fasta },
                STRUCTURAL_ANNOTATION.out.gbk,
                ch_submission_template,
                STRUCTURAL_ANNOTATION.out.busco_genome_summary.map { meta, summary -> summary },
                params.species,
                params.strain,
                true,
                true
            )
        } else {
            log.info "Skipping functional annotation (--no_functional_anno set)."
        }
    }

}

// ── SETUP entry point: download all pipeline databases ────────────────────────
workflow SETUP {
    DOWNLOAD_DATABASES()
}

// ── Pipeline completion hook ──────────────────────────────────────────────────
workflow.onComplete {
    log.info ""
    if (params.db_dir) {
        if (workflow.success) {
            log.info "Databases downloaded to: ${params.db_dir}"
            log.info ""
            log.info "Use these flags in your pipeline run:"
            log.info "  --databases           ${params.db_dir}"
            log.info "  --dbcan_db            ${params.db_dir}/dbcan"
            log.info "  --merops_db           ${params.db_dir}/merops/merops_pepunit.dmnd"
            log.info "  --phi_base_db         ${params.db_dir}/phi_base/phi_base.dmnd"
            log.info "  --IPS6_databases_path ${params.db_dir}/interproscan"
        } else {
            log.error "Database setup failed. See above for error details."
        }
    } else {
        if (workflow.success) {
            log.info "bagRNA completed successfully!"
            log.info "Results published to: ${params.outdir}"
        } else {
            log.error "bagRNA pipeline failed. See above for error details."
        }
    }
    log.info "Duration : ${workflow.duration}"
    log.info "CPU hours: ${workflow.stats.computeTimeFmt ?: '(not available)'}"
    log.info ""
}
