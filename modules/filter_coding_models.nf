process FILTER_CODING_MODELS {
    tag "filter_coding"
    label 'process_low'
    container 'ubuntu:22.04'

    publishDir "${params.outdir}/mikado", mode: 'copy'

    input:
    path mikado_pick_gff

    output:
    path "coding_mikado_models.gff", emit: coding_gff

    script:
    """
    # Identify ncRNA and mRNA without start/stop codons — collect their parent IDs.
    # Parent= is matched by prefix, not by fixed field position, since attribute
    # order in column 9 is not guaranteed across Mikado versions.
    awk 'BEGIN{FS=OFS="\\t"}{
        if((\$3=="ncRNA") ||
           (\$3=="mRNA" && (\$9~"has_start_codon=False" || \$9~"has_stop_codon=False"))) {
            n=split(\$9, attrs, ";");
            for (i=1; i<=n; i++) { if (attrs[i] ~ /^Parent=/) { print substr(attrs[i], 8); break } }
        }
    }' ${mikado_pick_gff} \\
        > bad_model_ids.txt

    # Filter those IDs out of the GFF and remove superlocus lines
    grep -w -F -v -f bad_model_ids.txt ${mikado_pick_gff} \\
        | sed 's/Name=.*//g' \\
        | awk 'BEGIN{FS=OFS="\\t"}{
            if(\$3!="superlocus") print
          }' \\
        > coding_mikado_models.gff
    """

    stub:
    """
    touch coding_mikado_models.gff
    """
}
