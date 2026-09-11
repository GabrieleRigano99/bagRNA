process GENERATE_REPORT {
    tag "generate_report"
    label 'process_low'
    // Chromium baked in at build time (docker/generate_report/Dockerfile) —
    // installing it at container RUNTIME (previously: `apt-get install
    // chromium` in this script) silently never worked in the real pipeline:
    // Nextflow always runs containers as the host's non-root UID, so
    // apt-get always failed with "Permission denied" (best-effort PDF
    // render meant this failed silently — confirmed missing on both the
    // Cryptococcus depauperatus and Syncephalastrum racemosum runs,
    // 2026-09-09, before this fix).
    container 'gabrielerigano/bagrna-report:latest'
    containerOptions '-v /usr/bin/ps:/usr/bin/ps'  // base image has no ps; nextflow needs it for task metrics

    publishDir "${params.outdir}/annotation_report", mode: 'copy'

    input:
    path genome_fasta
    path final_gff,          stageAs: 'input_final.gff'
    path proteins_faa
    path busco_summary,      stageAs: 'input_busco_summary.txt'
    path annotation_stats,   stageAs: 'input_annotation_stats.txt'
    path dbcan_overview,     stageAs: 'input_dbcan_overview.tsv'
    path antismash_json,     stageAs: 'input_antismash.json'
    path phi_base_tsv,       stageAs: 'input_phi_base.tsv'
    path eggnog_annotations, stageAs: 'input_eggnog.annotations'
    path kegg_module_json,   stageAs: 'input_kegg_modules.json'
    path cog_def_tab
    val  species
    val  strain

    output:
    path "*_annotation_summary.tbl",  emit: tbl
    path "*_annotation_report.html", emit: html
    path "*_annotation_report.pdf",  emit: pdf, optional: true

    script:
    // Every input below besides genome_fasta/proteins_faa/cog_def_tab may be
    // the assets/NO_FILE sentinel (that pipeline stage was skipped) — .name
    // still reflects the ORIGINAL (pre-stageAs) filename here, matching the
    // same check already used in modules/kegg_ko_info.nf.
    def opt = { flag, f, staged -> f.name != 'NO_FILE' ? "--${flag} ${staged}" : '' }
    """
    generate_annotation_report.py \\
        --genome-fasta ${genome_fasta} \\
        ${opt('final-gff', final_gff, 'input_final.gff')} \\
        --proteins-faa ${proteins_faa} \\
        ${opt('busco-summary', busco_summary, 'input_busco_summary.txt')} \\
        ${opt('annotation-stats', annotation_stats, 'input_annotation_stats.txt')} \\
        ${opt('dbcan-overview', dbcan_overview, 'input_dbcan_overview.tsv')} \\
        ${opt('antismash-json', antismash_json, 'input_antismash.json')} \\
        ${opt('phi-base-tsv', phi_base_tsv, 'input_phi_base.tsv')} \\
        ${opt('eggnog-annotations', eggnog_annotations, 'input_eggnog.annotations')} \\
        ${opt('kegg-module-completeness', kegg_module_json, 'input_kegg_modules.json')} \\
        --cog-def-tab ${cog_def_tab} \\
        --species ${species} \\
        --strain ${strain} \\
        --outdir .
    """

    stub:
    """
    touch ${species}_${strain}_annotation_summary.tbl
    touch ${species}_${strain}_annotation_report.html
    touch ${species}_${strain}_annotation_report.pdf
    """
}
