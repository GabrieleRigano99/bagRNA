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
    # AGAT's --type does exact string matching against column 3, not SO
    # ontology expansion -- "ncRNA" never matches the literal "tRNA"/"rRNA"
    # feature types this pipeline's GFFs actually use (from Barrnap/tRNAscan),
    # so a single --type ncRNA call always produced an empty file. Extract
    # each real type separately and concatenate.
    touch /data/ncrna_transcripts.fasta
    for t in tRNA rRNA; do
        agat_sp_extract_sequences.pl \\
            --gff /data/${final_gff} \\
            --fasta /data/${fasta} \\
            --type \$t \\
            --merge \\
            -o /data/ncrna_\${t}.fasta || true
        [ -s /data/ncrna_\${t}.fasta ] && cat /data/ncrna_\${t}.fasta >> /data/ncrna_transcripts.fasta
    done
    """

    stub:
    """
    touch ncrna_transcripts.fasta
    """
}
