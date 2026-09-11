process AGAT_FIX_CDS_LIFTOVER {
    tag "fix_liftover"
    label 'process_medium'
    containerOptions '-v $PWD:/data'
    container 'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1'

    publishDir "${params.outdir}/liftover", mode: 'copy'

    input:
    path lifted_gff
    path fasta

    output:
    path "liftover.gff", emit: liftover_gff

    script:
    """
    # Fix CDS phases of the lifted annotation
    agat_sp_fix_cds_phases.pl \\
        --gff /data/${lifted_gff} \\
        --fasta /data/${fasta} \\
        -o /data/fix_cds_lifted_anno.gff

    # Extract only records that have an associated mRNA (remove orphan features)
    awk 'BEGIN{FS=OFS="\\t"}{if(\$3~"mRNA") print \$9}' /data/fix_cds_lifted_anno.gff \\
        | sed 's/;/\\t/g' \\
        | sed 's/Parent=//g' \\
        | cut -f2 \\
        | grep -w -F -f - /data/fix_cds_lifted_anno.gff \\
        | sed 's/,""//g' \\
        | awk 'BEGIN{print "##gff-version 3"; FS=OFS="\\t"} {print}' \\
        > /data/liftover.gff
    """

    stub:
    """
    touch liftover.gff
    """
}
