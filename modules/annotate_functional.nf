// Merge all functional-annotation tool outputs into the final structural GFF3.
// Delegates to bin/merge_functional_annotations.py (single source of truth,
// also usable standalone) rather than an embedded heredoc, so the pipeline
// and manual runs share the exact same, tested merging logic.

process ANNOTATE_FUNCTIONAL {
    tag "annotate_functional"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.11'

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path gff,                stageAs: 'input.gff'
    path eggnog_anno,        stageAs: 'input_eggnog'
    path interpro_tsv,       stageAs: 'input_interpro'
    path kofamscan_tsv,      stageAs: 'input_kofamscan'
    path infernal_tbl,       stageAs: 'input_infernal'
    path effectorp3_txt,     stageAs: 'input_effectorp3'
    path dbcan_overview,     stageAs: 'input_dbcan_overview'
    path dbcan_cgc,          stageAs: 'input_dbcan_cgc'
    path dbcan_tc,           stageAs: 'input_dbcan_tc'
    path merops_tsv,         stageAs: 'input_merops'
    path phi_base_tsv,       stageAs: 'input_phi_base'
    path ko_info,            stageAs: 'input_ko_info'
    path go_obo,             stageAs: 'input_go_obo'
    path antismash_gbk,      stageAs: 'input_antismash'
    path gene2product_file,  stageAs: 'input_gene2product'
    path id_map,             stageAs: 'input_id_map'

    output:
    path "functional_annotation.tsv", emit: annotation_table
    path "annotated.gff3",            emit: annotated_gff
    path "annotation_stats.txt",      emit: stats

    script:
    """
    python3 ${workflow.projectDir}/bin/merge_functional_annotations.py \\
        --gff ${gff} \\
        --eggnog ${eggnog_anno} \\
        --interpro ${interpro_tsv} \\
        --kofamscan ${kofamscan_tsv} \\
        --infernal ${infernal_tbl} \\
        --effectorp3 ${effectorp3_txt} \\
        --dbcan-overview ${dbcan_overview} \\
        --dbcan-cgc ${dbcan_cgc} \\
        --dbcan-tc ${dbcan_tc} \\
        --merops ${merops_tsv} \\
        --phi-base ${phi_base_tsv} \\
        --ko-info ${ko_info} \\
        --go-obo ${go_obo} \\
        --antismash-gbk ${antismash_gbk} \\
        --gene2product ${gene2product_file} \\
        --id-map ${id_map} \\
        --out-gff annotated.gff3 \\
        --out-tsv functional_annotation.tsv \\
        --out-stats annotation_stats.txt
    """

    stub:
    """
    touch functional_annotation.tsv annotated.gff3 annotation_stats.txt
    """
}
