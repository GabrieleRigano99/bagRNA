process KOFAMSCAN {
    tag "kofamscan"
    label 'process_high'
    container 'quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2'

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path proteins_faa
    path profiles_dir
    path ko_list

    output:
    path "kofamscan_result.tsv", emit: kofamscan_tsv

    script:
    """
    exec_annotation \\
        -f detail-tsv \\
        --cpu ${task.cpus} \\
        --profile ${profiles_dir} \\
        --ko-list ${ko_list} \\
        -o raw_kofamscan.tsv \\
        ${proteins_faa}

    # Keep only significant hits (marked with * in column 1)
    awk '\$2 != "" && \$1 == "*"' raw_kofamscan.tsv > kofamscan_result.tsv
    """

    stub:
    """
    touch kofamscan_result.tsv
    """
}
