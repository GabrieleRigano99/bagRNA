process FETCH_KEGG_MODULE_COMPLETENESS {
    tag "kegg_module_completeness"
    label 'process_low'
    container 'python:3.11-slim'
    containerOptions '-v /usr/bin/ps:/usr/bin/ps'  // python:3.11-slim has no ps; nextflow needs it for task metrics

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path kegg_annotations_tsv
    path eggnog_annotations, stageAs: 'input_eggnog.annotations'

    output:
    path "kegg_module_completeness.json", emit: json

    script:
    // eggNOG often calls genes the KO pipeline misses; folding its KO/EC
    // calls in as extra evidence avoids understating real pathway
    // completeness (see bin/fetch_kegg_module_completeness.py's docstring).
    def eggnog_arg = eggnog_annotations.name != 'NO_FILE' ? "--eggnog input_eggnog.annotations" : ''
    """
    fetch_kegg_module_completeness.py ${kegg_annotations_tsv} kegg_module_completeness.json ${eggnog_arg}
    """

    stub:
    """
    echo '[]' > kegg_module_completeness.json
    """
}
