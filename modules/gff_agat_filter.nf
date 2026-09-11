process GFF_AGAT_FILTER {
    tag "gff_agat_filter"
    label 'process_high'
    container 'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1'

    input:
    path merged_gff
    path fasta
    val  max_gene_length

    output:
    path "longest_iso.gff", emit: filtered_gff

    script:
    """
    agat_sp_fix_overlaping_genes.pl \\
        --gff ${merged_gff} \\
        -o fix_olp.gff

    agat_sp_fix_features_locations_duplicated.pl \\
        --gff fix_olp.gff \\
        -o fix_locdup.gff

    agat_sp_filter_incomplete_gene_coding_models.pl \\
        -f ${fasta} \\
        --gff fix_locdup.gff \\
        -o filter_incomplete_gene.gff

    agat_sp_filter_by_ORF_size.pl \\
        -s 49 \\
        --gff filter_incomplete_gene.gff \\
        -o orf_filt_incomplete.gff

    # The container's /bin/awk is BusyBox awk, not gawk: its match() silently
    # ignores a 3rd (capture-array) argument instead of erroring, so the old
    # match-with-capture-array form left the ID variable permanently empty.
    # That blank ID broke this line's field count, so the length field in the
    # second awk was always empty and the length comparison never matched —
    # kill_list.txt has been silently empty on every run in this pipeline's
    # history, meaning --max_gene_length has never actually filtered anything.
    # Extract the ID via split()/substr() instead, which BusyBox awk supports.
    awk '\$3=="gene"{
        n=split(\$9, attrs, ";");
        id="";
        for (i=1; i<=n; i++) { if (attrs[i] ~ /^ID=/) { id=substr(attrs[i], 4); break } }
        print id, \$5-\$4+1
    }' orf_filt_incomplete_sup49.gff \\
        | sort -k2,2nr \\
        | awk -v max="${max_gene_length}" '\$2 >= max {print}' \\
        | cut -f1 -d " " \\
        > kill_list.txt

    agat_sp_filter_feature_from_kill_list.pl \\
        --gff orf_filt_incomplete_sup49.gff \\
        --kill_list kill_list.txt \\
        -o longest_iso.gff
    """

    stub:
    """
    touch longest_iso.gff
    """
}
