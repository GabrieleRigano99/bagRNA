// functional_annotation.nf
// Functional annotation sub-workflow: stages 41-50
// Covers InterProScan 6 through the ANNOTATE_FUNCTIONAL merge

nextflow.enable.dsl = 2

include { PREPARE_INTERPROSCAN;
          INTERPROSCAN as INTERPROSCAN6 } from '../subworkflows/interproscan6/workflows/interproscan'
include { EGGNOG            } from '../modules/eggnog'
include { DEEPKOALA         } from '../modules/deepkoala'
include { KEGG_ANNOTATE     } from '../modules/kegg_annotate'
include { KEGG_KO_INFO      } from '../modules/kegg_ko_info'
include { INFERNAL          } from '../modules/infernal'
include { SIGNALP6          } from '../modules/signalp6'
include { PHOBIUS           } from '../modules/phobius'
include { EFFECTORP3        } from '../modules/effectorp3'
include { ANTISMASH         } from '../modules/antismash'
include { MEROPS                  } from '../modules/merops'
include { PHI_BASE                } from '../modules/phi_base'
include { RUN_DBCAN               } from '../modules/run_dbcan'
include { ANNOTATE_FUNCTIONAL     } from '../modules/annotate_functional'
include { GFF3_TO_TBL;
          NCBI_SUBMISSION         } from '../modules/ncbi_submission'
include { FETCH_KEGG_MODULE_COMPLETENESS } from '../modules/fetch_kegg_module_completeness'
include { GENERATE_REPORT                } from '../modules/generate_report'

// Resolve a database path: if it's a directory, return the .dmnd file inside
// (or FASTA if no .dmnd). Ensures Nextflow stages the real file so Docker
// auto-mounts the containing path rather than a broken symlink chain.
def resolveDbFile(p) {
    if (!p.isDirectory()) return p
    def dmnd  = p.listFiles()?.find { it.name.endsWith('.dmnd') }
    def fasta = p.listFiles()?.find { it.name =~ /\.(fasta|fas|fa)$/ }
    return dmnd ?: fasta ?: p
}

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_proteins_faa        // path: final protein FASTA
    ch_ncrna_fasta         // path: ncRNA transcript FASTA
    ch_final_gff           // path: final structural GFF (or sentinel)
    ch_fasta               // path: genome FASTA
    ch_gbk                 // path: table2asn .gbk files
    ch_submission_template // path: NCBI .sbt template
    ch_busco_summary       // path: compleasm genome-mode summary.txt (or NO_FILE sentinel)
    species                // val
    strain                 // val
    run_funannotate        // val(boolean): false skips ANNOTATE_FUNCTIONAL (no --final_gff provided)
    run_infernal           // val(boolean): false skips INFERNAL (no ncRNA FASTA provided)
    gbk_provided           // val(boolean): false skips ANTISMASH (no real .gbk, ch_gbk is the NO_FILE sentinel)

    main:

    def no_file = file("${projectDir}/assets/NO_FILE")

    // Resolve database paths from params.databases
    def db_dir       = params.databases ? file(params.databases) : null
    def rfam_cm      = db_dir ? file("${params.databases}/Rfam.cm")     : null
    def rfam_clanin  = db_dir ? file("${params.databases}/Rfam.clanin") : null
    // Optional databases: explicit flag takes priority, then fall back to --databases subdir
    def dbcan_path    = params.dbcan_db    ? file(params.dbcan_db)
                      : db_dir            ? file("${params.databases}/dbcan")
                      : null
    def merops_path   = params.merops_db   ? file(params.merops_db)
                      : db_dir            ? resolveDbFile(file("${params.databases}/merops"))
                      : null
    def phi_base_path = params.phi_base_db ? file(params.phi_base_db)
                      : db_dir            ? resolveDbFile(file("${params.databases}/phi_base"))
                      : null
    def go_obo_path        = params.go_obo       ? file(params.go_obo)
                            : db_dir             ? file("${params.databases}/go-basic.obo")
                            : null
    def go_obo_file       = (go_obo_path && go_obo_path.exists()) ? go_obo_path : no_file
    def gene2product_path  = params.gene2product ? file(params.gene2product)
                            : db_dir             ? file("${params.databases}/ncbi_cleaned_gene_products.txt")
                            : null
    def gene2product_file = (gene2product_path && gene2product_path.exists()) ? gene2product_path : no_file
    def id_map_file       = params.id_map       ? file(params.id_map)       : no_file
    def eggnog_data_dir    = params.eggnog_data_dir ? params.eggnog_data_dir
                            : db_dir                ? "${params.databases}/eggnog_data"
                            : null
    def eggnog_data_dir_available = eggnog_data_dir && file(eggnog_data_dir).exists()

    // ── 41. InterProScan 6 (uses InterPro API by default; set --IPS6_databases_path for local mode) ──
    if (!params.no_interpro) {
        def ips6_apps_config = "${projectDir}/subworkflows/interproscan6/conf/applications.config"
        def ips6_apps        = params.interproscan6_apps ? params.interproscan6_apps.tokenize(',') : []
        def ips6_fallback_dir = db_dir ? file("${params.databases}/interproscan") : null
        def ips6_data_dir    = params.IPS6_databases_path
                                   ? file(params.IPS6_databases_path).toAbsolutePath().toString()
                                   : (ips6_fallback_dir && ips6_fallback_dir.exists())
                                   ? ips6_fallback_dir.toAbsolutePath().toString()
                                   : null
        def ips6_outprefix   = "${params.outdir}/functional_annotation/interpro_output"
        file("${params.outdir}/functional_annotation").mkdirs()

        PREPARE_INTERPROSCAN(ips6_apps_config, ips6_apps, params.use_gpu)
        def ips6_config = PREPARE_INTERPROSCAN.out.apps_config.val

        // Link local Phobius installation
        if (params.phobius_path) {
            ips6_config['phobius']['dir'] = params.phobius_path
            if (ips6_apps && !ips6_apps.contains('phobius'))
                ips6_apps = ips6_apps + ['phobius']
        }

        // Link local SignalP installation (eukarya mode)
        if (params.signalp_path && !params.no_signalp) {
            ips6_config['signalp_euk']['dir'] = params.signalp_path
            if (ips6_apps && !ips6_apps.contains('signalp_euk'))
                ips6_apps = ips6_apps + ['signalp_euk']
        }

        // [] means "no apps" in MatchesApiClient — derive full list from config instead
        if (!ips6_apps) {
            if (!ips6_data_dir) {
                // No local databases: limit to apps that can run without db files.
                // This makes matches_api_apps empty in PREPARE_APPLICATIONS, which
                // skips SCAN_REMAINING (avoiding null-dirpath failures for API misses).
                ips6_apps = ['coils', 'mobidblite', 'tmbed']
                if (params.phobius_path) ips6_apps += ['phobius']
                if (params.signalp_path && !params.no_signalp) ips6_apps += ['signalp_euk']
                log.warn "InterProScan6: no '--IPS6_databases_path' provided — limiting to local-only apps: ${ips6_apps.join(', ')}. Provide --IPS6_databases_path for full analysis."
            } else {
                def (allApps, appsErr, appsWarn) = uk.ac.ebi.interpro.InterProScan.validateApplications(
                    null, null, ips6_config, true
                )
                if (!allApps) {
                    log.error "InterProScan app validation failed: ${appsErr}"
                    System.exit(1)
                }
                if (appsWarn) log.warn "InterProScan: ${appsWarn}"
                ips6_apps = allApps
            }
        }

        // Re-apply GPU flag after apps list is finalised — PREPARE_INTERPROSCAN ran
        // enableGpuAcceleration with an empty list (before auto-selection), so tmbed
        // and signalp_euk still had use_gpu=false in the returned config.
        if (params.use_gpu) {
            def (gpuConfig, gpuWarn) = uk.ac.ebi.interpro.InterProScan.enableGpuAcceleration(
                true, ips6_apps, ips6_config
            )
            ips6_config = gpuConfig
            if (gpuWarn) log.warn "InterProScan: ${gpuWarn}"
            // TMBed GPU batch = batch_size * 10; reduce to avoid GPU OOM.
            // Default 3000 * 10 = 30000 exceeds most GPU VRAM budgets.
            if (ips6_config.containsKey('tmbed') && ips6_config.tmbed.use_gpu) {
                ips6_config['tmbed']['batch_size'] = params.interproscan6_tmbed_signalp6_gpu_batch_size ?: 200
            }
        }

        INTERPROSCAN6(
            ch_proteins_faa,
            ips6_apps,
            ips6_config,
            ips6_data_dir,
            ips6_outprefix,
            ["xml", "tsv"],
            params.interproscan6_interpro_version?.toString() ?: "latest",
            "6.0.0",                // interproscan_version
            "InterProScan6",        // interproscan_name
            params.interproscan6_no_matches_api ?: false,
            "https://www.ebi.ac.uk/interpro/matches/api",
            100,                    // matches_api_chunk_size
            3,                      // matches_api_max_retries
            5000,                   // batch_size
            500,                    // sub_batch_size
            false,                  // nucleic
            false,                  // skip_repr_locations
            params.interproscan6_goterms  ?: false,
            params.interproscan6_pathways ?: false,
            false,                  // globus
            false                   // enforce_compatibility
        )

        ch_interpro_xml = INTERPROSCAN6.out.output_files
            .flatten()
            .filter { it.toString().endsWith('.xml') }
            .map    { file(it.toString()) }
        ch_interpro_tsv = INTERPROSCAN6.out.output_files
            .flatten()
            .filter { it.toString().endsWith('.tsv') }
            .map    { file(it.toString()) }
    } else {
        ch_interpro_xml = Channel.value(no_file)
        ch_interpro_tsv = Channel.value(no_file)
    }

    // eggNOG-mapper v7 data dir: explicit flag first, else --databases/eggnog_data
    if (!params.no_eggnog && eggnog_data_dir_available) {
        EGGNOG(ch_proteins_faa, species, strain, file(eggnog_data_dir))
        ch_eggnog_annotations = EGGNOG.out.eggnog_annotations
    } else {
        if (!params.no_eggnog) {
            log.warn "Skipping EGGNOG: --eggnog_data_dir is required (or set --no_eggnog to silence this)."
        }
        ch_eggnog_annotations = Channel.value(no_file)
    }

    // ── 42–44. Database-dependent analyses (skipped when --databases absent) ──
    if (params.databases) {
        // KO assignment via DeepKOALA (GRU deep-learning classifier) in place
        // of KofamScan's HMM search — benchmarked faster (60-107x on GPU,
        // ~3.5x even CPU-only) and higher-recall at comparable precision
        // (cross-validated against eggNOG's independent KO calls on both
        // Candida orthopsilosis and Sporothrix schenckii). Output is
        // reformatted to KofamScan's own TSV layout so every downstream
        // consumer is unchanged; see modules/deepkoala.nf for the field
        // mapping. modules/kofamscan.nf is kept on disk, unused, in case of
        // a future revert.
        DEEPKOALA(ch_proteins_faa)
        ch_kofamscan_tsv = DEEPKOALA.out.kofamscan_tsv

        KEGG_ANNOTATE(ch_kofamscan_tsv)
        ch_kegg_tsv = KEGG_ANNOTATE.out.kegg_tsv

        // KO-keyed KEGG lookup covering KofamScan + eggNOG KOs (for gene
        // symbols / KO products in ANNOTATE_FUNCTIONAL), cached off KEGG_ANNOTATE.
        KEGG_KO_INFO(ch_kofamscan_tsv, ch_eggnog_annotations, ch_kegg_tsv)
        ch_ko_info = KEGG_KO_INFO.out.ko_info

        if (run_infernal) {
            INFERNAL(ch_ncrna_fasta, rfam_cm, rfam_clanin)
            ch_infernal_table = INFERNAL.out.infernal_table
        } else {
            ch_infernal_table = Channel.value(no_file)
        }
    } else {
        ch_kofamscan_tsv      = Channel.value(no_file)
        ch_kegg_tsv           = Channel.value(no_file)
        ch_infernal_table     = Channel.value(no_file)
    }

    // ── 45. SignalP6 signal peptide — standalone only when interproscan6 disabled ──
    if (!params.no_signalp && params.signalp_path && params.no_interpro) {
        SIGNALP6(ch_proteins_faa)
        ch_signalp_output = SIGNALP6.out.signalp_results
    } else {
        if (!params.no_signalp && !params.signalp_path) {
            log.warn "Skipping SignalP6: '--signalp_path' not provided."
        }
        ch_signalp_output = Channel.value(no_file)
    }

    // ── 46. Phobius — standalone only when interproscan6 disabled ────────────
    if (params.phobius_path && params.no_interpro) {
        PHOBIUS(ch_proteins_faa)
        ch_phobius_output = PHOBIUS.out.phobius_results
    } else {
        ch_phobius_output = Channel.value(no_file)
    }

    // ── 47. EffectorP-3 effector prediction (optional) ───────────────────────
    if (!params.no_effectorp3) {
        EFFECTORP3(ch_proteins_faa)
        ch_effectorp3_output = EFFECTORP3.out.effectorp3_results
    } else {
        ch_effectorp3_output = Channel.value(no_file)
    }

    // ── 48a. run_dbcan CAZyme + CGC/PUL annotation (optional) ────────────────
    if (dbcan_path && !params.no_dbcan) {
        RUN_DBCAN(ch_proteins_faa, ch_final_gff, dbcan_path)
        ch_dbcan_overview = RUN_DBCAN.out.overview
        ch_dbcan_cgc      = RUN_DBCAN.out.cgc_standard
        ch_dbcan_gff      = RUN_DBCAN.out.cgc_gff
    } else {
        ch_dbcan_overview = Channel.value(no_file)
        ch_dbcan_cgc      = Channel.value(no_file)
        ch_dbcan_gff      = Channel.value(no_file)
    }
    // TCDB transporter hits (diamond.out.tc) sit alongside overview.tsv in the
    // dbcan output dir (retained by RUN_DBCAN's `dbcan_out/*` output) — reach
    // it as a sibling rather than re-running RUN_DBCAN just to add an emit.
    ch_dbcan_tc = ch_dbcan_overview.map { f ->
        def tc = file("${f.parent}/diamond.out.tc")
        tc.exists() ? tc : no_file
    }

    // ── 48b. MEROPS peptidase annotation (optional) ──────────────────────────
    if (merops_path && !params.no_merops) {
        MEROPS(ch_proteins_faa, merops_path)
        ch_merops_tsv = MEROPS.out.merops_tsv
    } else {
        ch_merops_tsv = Channel.value(no_file)
    }

    // ── 48b. PHI-base pathogen-host interaction (optional) ───────────────────
    if (phi_base_path && !params.no_phi_base) {
        PHI_BASE(ch_proteins_faa, phi_base_path)
        ch_phi_base_tsv = PHI_BASE.out.phi_base_tsv
    } else {
        ch_phi_base_tsv = Channel.value(no_file)
    }

    // ── 49. AntiSMASH secondary metabolite clusters (optional) ───────────────
    if (!params.no_antismash && gbk_provided) {
        ANTISMASH(ch_gbk.first(), species, strain)
        ch_antismash_gbk = ANTISMASH.out.antismash_dir
            .map { dir -> file("${dir}/antismash_${species}_${strain}.gbk") }
        ch_antismash_json = ANTISMASH.out.antismash_dir
            .map { dir -> file("${dir}/antismash_${species}_${strain}.json") }
    } else {
        ch_antismash_gbk  = Channel.value(no_file)
        ch_antismash_json = Channel.value(no_file)
    }

    // ── 50a. bagRNA annotation merger: all sources → TSV + annotated GFF3 ───────
    if (run_funannotate) {
        ANNOTATE_FUNCTIONAL(
            ch_final_gff,
            ch_eggnog_annotations,
            ch_interpro_tsv,
            ch_kofamscan_tsv,
            ch_infernal_table,
            ch_effectorp3_output,
            ch_dbcan_overview,
            ch_dbcan_cgc,
            ch_dbcan_tc,
            ch_merops_tsv,
            ch_phi_base_tsv,
            ch_ko_info,
            go_obo_file,
            ch_antismash_gbk,
            gene2product_file,
            id_map_file
        )
        ch_annotation_table = ANNOTATE_FUNCTIONAL.out.annotation_table
        ch_annotated_gff    = ANNOTATE_FUNCTIONAL.out.annotated_gff
        ch_annotation_stats = ANNOTATE_FUNCTIONAL.out.stats

        // ── 50b. NCBI submission files (gff3_to_tbl.py + table2asn) ──────────
        if (!params.no_ncbi_submission) {
            GFF3_TO_TBL(
                ch_annotated_gff,
                ch_annotation_table,
                species,
                strain
            )
            NCBI_SUBMISSION(
                ch_fasta,
                GFF3_TO_TBL.out.tbl,
                ch_submission_template,
                species,
                strain,
                params.codon_table,
                params.locus_tag
            )
        }
    } else {
        ch_annotation_table = Channel.value(no_file)
        ch_annotated_gff    = Channel.value(no_file)
        ch_annotation_stats = Channel.value(no_file)
    }

    // ── 51. KEGG module completeness — live KEGG REST lookups keyed on the
    //        modules already present in KEGG_ANNOTATE's output, so this can
    //        only run when that ran (i.e. --databases was given). ──────────
    if (params.databases && !params.no_report) {
        FETCH_KEGG_MODULE_COMPLETENESS(ch_kegg_tsv, ch_eggnog_annotations)
        ch_kegg_module_completeness = FETCH_KEGG_MODULE_COMPLETENESS.out.json
    } else {
        ch_kegg_module_completeness = Channel.value(no_file)
    }

    // ── 52. Annotation report: PDF + summary table over every stage above.
    //        Every input besides the genome/proteins/COG-reference asset is
    //        optional — sections for skipped stages are simply left out
    //        rather than failing the whole report (see generate_report.nf).
    if (!params.no_report) {
        def cog_def_tab = file("${projectDir}/assets/cog-20.def.tab")
        GENERATE_REPORT(
            ch_fasta,
            ch_final_gff,
            ch_proteins_faa,
            ch_busco_summary,
            ch_annotation_stats,
            ch_dbcan_overview,
            ch_antismash_json,
            ch_phi_base_tsv,
            ch_eggnog_annotations,
            ch_kegg_module_completeness,
            cog_def_tab,
            species,
            strain
        )
    }


    emit:
    interpro_xml        = ch_interpro_xml
    interpro_tsv        = ch_interpro_tsv
    eggnog_annotations  = ch_eggnog_annotations
    kofamscan_tsv       = ch_kofamscan_tsv
    kegg_tsv            = ch_kegg_tsv
    infernal_table      = ch_infernal_table
    signalp_results     = ch_signalp_output
    effectorp3_results  = ch_effectorp3_output
    merops_tsv          = ch_merops_tsv
    phi_base_tsv        = ch_phi_base_tsv
    dbcan_overview      = ch_dbcan_overview
    dbcan_cgc           = ch_dbcan_cgc
    dbcan_gff           = ch_dbcan_gff
    annotation_table    = ch_annotation_table
    annotated_gff       = ch_annotated_gff
    annotation_stats    = ch_annotation_stats
}
