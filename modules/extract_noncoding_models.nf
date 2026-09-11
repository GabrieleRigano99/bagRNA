process EXTRACT_NONCODING_MODELS {
    tag "extract_noncoding"
    label 'process_low'
    container 'ubuntu:22.04'

    publishDir "${params.outdir}/mikado", mode: 'copy'

    input:
    path mikado_pick_gff

    output:
    path "noncoding_mikado_models.gff", emit: noncoding_gff

    script:
    """
    # Collect IDs of ncRNA features and mRNA lacking start/stop codons
    awk 'BEGIN{FS=OFS="\\t"}{
        if((\$3=="ncRNA") ||
           (\$3=="mRNA" && (\$9~"has_start_codon=False" || \$9~"has_stop_codon=False")))
            print \$9
    }' ${mikado_pick_gff} \\
        | sed 's/;/\\t/g' \\
        | sed 's/Parent=//g' \\
        | cut -f2 \\
        > noncoding_ids.txt

    # Extract those records from the GFF and reclassify mRNA → ncRNA
    grep -w -F -f noncoding_ids.txt ${mikado_pick_gff} \\
        | awk 'BEGIN{print "##gff-version 3"; FS=OFS="\\t"}{
            if(\$9!~"has_start_codon=True" || \$9!~"has_stop_codon=True") print
          }' \\
        | awk 'BEGIN{FS=OFS="\\t"}{if(\$3!~"UTR" && \$3!~"CDS") print}' \\
        | sed 's/Name=.*//g' \\
        | sed 's/mRNA/ncRNA/g' \\
        | awk 'BEGIN{FS=OFS="\\t"}{if(\$3!~"UTR" && \$3!~"CDS") print}' \\
        | sed 's/alias.*//g' \\
        | sed 's/ncRNA_gene/gene/g' \\
        | sed 's/fpkm=.*//g' \\
        > noncoding_mikado_models.gff
    """

    stub:
    """
    touch noncoding_mikado_models.gff
    """
}
