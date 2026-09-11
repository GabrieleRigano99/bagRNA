process INFERNAL {
    tag "infernal"
    label 'process_high'
    container 'quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_2'

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path ncrna_fasta
    path rfam_cm
    path rfam_clanin

    output:
    path "infernal_table.txt", emit: infernal_table

    script:
    """
    cmpress ${rfam_cm}

    cmscan \\
        --cut_ga \\
        --rfam \\
        --nohmmonly \\
        --fmt 2 \\
        --tblout infernal_table.txt \\
        --cpu ${task.cpus} \\
        --oskip \\
        --clanin ${rfam_clanin} \\
        ${rfam_cm} \\
        ${ncrna_fasta}
    """

    stub:
    """
    touch infernal_table.txt
    """
}
