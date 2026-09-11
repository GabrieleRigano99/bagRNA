process GFF_CLEAN_FILTER {
    tag "gff_clean_filter"
    label 'process_high'
    container 'nextgenusfs/funannotate:v1.8.17'

    publishDir "${params.outdir}/structural_annotation", mode: 'copy'

    input:
    path renamed_gff   // already locus-tag-renamed (RENAME_LOCUS_TAG_IDS)
    path fasta
    val  species
    val  strain
    val  locus_tag

    output:
    path "genes_kept.gff", emit: filtered_gff

    script:
    """
    # Discrepancy removal via funannotate gff2tbl / tbl2gbk
    funannotate util gff2tbl \\
        -g ${renamed_gff} \\
        -f ${fasta} \\
        > tmp_sorted_structural.tbl

    funannotate util tbl2gbk \\
        -i tmp_sorted_structural.tbl \\
        -f ${fasta} \\
        -s "${species}" \\
        -o tmp_sorted_structural \\
        > genes_to_fix_or_remove.txt 2>&1 || true

    # Remove problematic gene models identified by tbl2gbk
    tail -n +3 genes_to_fix_or_remove.txt \\
        | cut -f1 \\
        | grep -v -F -w -f - ${renamed_gff} \\
        > genes_kept.gff || cp ${renamed_gff} genes_kept.gff

    # Safety guard: the word-boundary grep above can catastrophically match a
    # locus-tag fragment and remove EVERY gene — e.g. an invalid hyphenated
    # prefix like "1099-18" makes funannotate report the bare word "1099",
    # which then matches all "1099-18_*" IDs. Abort loudly here rather than
    # emit an empty annotation that only fails cryptically downstream at
    # table2asn ("Unable to load annotations").
    n_in=\$(awk -F'\\t' '\$3=="gene"' ${renamed_gff} | wc -l)
    n_out=\$(awk -F'\\t' '\$3=="gene"' genes_kept.gff | wc -l)
    if [ "\$n_in" -gt 0 ] && [ "\$n_out" -eq 0 ]; then
        echo "ERROR: gene-model cleanup removed ALL \$n_in genes (0 kept)." >&2
        echo "Most likely the locus_tag prefix '${locus_tag}' is malformed:" >&2
        echo "NCBI requires 3-12 alphanumeric chars, starting with a letter, no hyphens." >&2
        exit 1
    fi
    """

    stub:
    """
    touch genes_kept.gff
    """
}
