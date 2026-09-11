process PORTCULLIS {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/portcullis:1.2.4--py36h6f1349a_0'

    publishDir "${params.outdir}/portcullis/${meta.id}", mode: 'copy', pattern: "Portcullis_${meta.id}/3-filt/*.bed"
    publishDir "${params.outdir}/portcullis/${meta.id}", mode: 'copy', pattern: "Portcullis_${meta.id}/portcullis.filtered.bam"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    val   orientation
    val   strandedness

    output:
    tuple val(meta), path("Portcullis_${meta.id}/3-filt/portcullis_filtered.pass.junctions.bed"), emit: junctions
    tuple val(meta), path("Portcullis_${meta.id}/portcullis.filtered.bam"),                       emit: filtered_bam

    script:
    """
    portcullis full \\
        --orientation ${orientation} \\
        --strandedness ${strandedness} \\
        --max_length ${params.max_intron_length} \\
        --min_cov 5 \\
        --bam_filter \\
        --output Portcullis_${meta.id} \\
        --threads ${task.cpus} \\
        ${fasta} \\
        ${bam}
    """

    stub:
    """
    mkdir -p Portcullis_${meta.id}/3-filt
    touch Portcullis_${meta.id}/3-filt/portcullis_filtered.pass.junctions.bed
    touch Portcullis_${meta.id}/portcullis.filtered.bam
    """
}
