process PHOBIUS {
    tag "phobius"
    label 'process_medium'
    container 'quay.io/biocontainers/perl:5.26.2'
    containerOptions "-v ${params.phobius_path}:/opt/phobius"

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path proteins_faa

    output:
    path "phobius_output.txt", emit: phobius_results

    script:
    """
    mkdir phobius_run
    cp /opt/phobius/phobius.pl /opt/phobius/phobius.options /opt/phobius/phobius.model phobius_run/
    cp /opt/phobius/decodeanhmm.64bit phobius_run/decodeanhmm
    chmod +x phobius_run/decodeanhmm

    perl phobius_run/phobius.pl \\
        -short \\
        ${proteins_faa} \\
        > phobius_output.txt
    """

    stub:
    """
    touch phobius_output.txt
    """
}
