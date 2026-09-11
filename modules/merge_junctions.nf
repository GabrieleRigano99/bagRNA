process MERGE_JUNCTIONS {
    tag "merge_junctions"
    label 'process_low'
    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1'

    publishDir "${params.outdir}/junctions", mode: 'copy'

    input:
    path junction_beds

    output:
    path "all_splice_junctions.bed", emit: merged_bed

    script:
    """
    cat ${junction_beds} \\
        | bedtools sort -i - \\
        | bedtools merge -i - \\
        > all_splice_junctions.bed
    """

    stub:
    """
    touch all_splice_junctions.bed
    """
}
