// Build a KO-keyed KEGG lookup (gene symbol, name, pathways) covering KOs
// from BOTH KofamScan and eggNOG, via bin/fetch_kegg_ko_info.py (KEGG REST).
// Reuses KEGG_ANNOTATE's kegg_annotations.tsv as a cache so already-fetched
// KOs are not re-requested. Feeds ANNOTATE_FUNCTIONAL's --ko-info.

process KEGG_KO_INFO {
    tag "kegg_ko_info"
    label 'process_low'
    container 'python:3.11-slim'
    containerOptions '-v /usr/bin/ps:/usr/bin/ps'

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path kofamscan_tsv, stageAs: 'kofamscan_result.tsv'
    path eggnog_anno,   stageAs: 'input_eggnog'
    path kegg_cache,    stageAs: 'kegg_annotations.tsv'

    output:
    path "ko_info.tsv", emit: ko_info

    script:
    def cache_arg = kegg_cache.name != 'NO_FILE' ? "--cache kegg_annotations.tsv" : ''
    """
    python3 ${workflow.projectDir}/bin/fetch_kegg_ko_info.py \\
        --kofamscan kofamscan_result.tsv \\
        --eggnog input_eggnog \\
        ${cache_arg} \\
        --out ko_info.tsv
    """

    stub:
    """
    touch ko_info.tsv
    """
}
