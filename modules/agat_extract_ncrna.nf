process AGAT_EXTRACT_NCRNA {
    tag "extract_ncrna"
    label 'process_medium'
    containerOptions '-v $PWD:/data'
    container 'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1'

    publishDir "${params.outdir}/structural_annotation", mode: 'copy'

    input:
    path final_gff
    path fasta

    output:
    path "ncrna_transcripts.fasta", emit: ncrna_fasta

    script:
    """
    agat_sp_extract_sequences.pl \\
        --gff /data/${final_gff} \\
        --fasta /data/${fasta} \\
        --type ncRNA \\
        --merge \\
        -o /data/ncrna_transcripts.fasta
    """

    stub:
    """
    touch ncrna_transcripts.fasta
    """
}
