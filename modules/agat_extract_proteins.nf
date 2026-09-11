process AGAT_EXTRACT_PROTEINS {
    tag "extract_proteins"
    label 'process_medium'
    containerOptions '-v $PWD:/data'
    container 'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1'

    publishDir "${params.outdir}/structural_annotation", mode: 'copy'

    input:
    path final_gff
    path fasta
    val  codon_table

    output:
    path "final_proteins.faa", emit: proteins_faa

    script:
    """
    agat_sp_extract_sequences.pl \\
        --gff /data/${final_gff} \\
        --fasta /data/${fasta} \\
        --aa \\
        --codon ${codon_table} \\
        --clean_internal_stop \\
        --clean_final_stop \\
        -o /data/final_proteins.faa
    """

    stub:
    """
    touch final_proteins.faa
    """
}
