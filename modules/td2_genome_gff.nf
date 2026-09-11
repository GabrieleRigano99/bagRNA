process TD2_GENOME_GFF {
    tag "td2_genome_gff"
    label 'process_low'
    container 'quay.io/biocontainers/transdecoder:5.7.0--pl5321hdfd78af_0'

    publishDir "${params.outdir}/transdecoder", mode: 'copy'

    input:
    path td2_gff3          // transcript-space GFF3 from TD2.Predict
    path mikado_gtf        // mikado_prepared.gtf (genome-space exon coordinates)
    path mikado_fasta      // mikado_prepared.fasta (transcript sequences)

    output:
    path "${td2_gff3.name.replaceAll(/\.gff3$/, '')}.genome.gff3", emit: genome_gff3

    script:
    def out = "${td2_gff3.name.replaceAll(/\.gff3$/, '')}.genome.gff3"
    """
    gtf_to_alignment_gff3.pl ${mikado_gtf} > alignment.gff3

    cdna_alignment_orf_to_genome_orf.pl \\
        ${td2_gff3} \\
        alignment.gff3 \\
        ${mikado_fasta} \\
        > ${out}
    """

    stub:
    def out = "${td2_gff3.name.replaceAll(/\.gff3$/, '')}.genome.gff3"
    """
    touch ${out}
    """
}
