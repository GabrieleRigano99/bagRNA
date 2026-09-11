process BUILD_EXTERNAL_SCORES {
    tag "external_scores"
    label 'process_low'
    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1'

    publishDir "${params.outdir}/mikado", mode: 'copy'

    input:
    path quant_sf
    path cpc2_result
    path busco_scores_tsv   // NO_FILE if absent

    output:
    path "external_scores_tpm_cpc2.tsv", emit: external_scores

    script:
    def has_busco = busco_scores_tsv.name != 'NO_FILE'
    """
    # Extract TPM column (col 4) from Salmon quant.sf, skip header, sort
    tail -n +2 ${quant_sf} \\
        | awk 'BEGIN{FS=OFS="\\t"}{print \$1, \$4}' \\
        | sort \\
        > tmp_tpm.tsv

    # Extract coding label column (col 7) from CPC2 result, skip header, sort
    tail -n +2 ${cpc2_result} \\
        | awk 'BEGIN{FS=OFS="\\t"}{print \$1, \$7}' \\
        | sort \\
        > tmp_cpc2.tsv

    # Add zero-TPM entries for transcripts present in CPC2 but absent from Salmon
    cut -f1 tmp_tpm.tsv \\
        | grep -v -F -w -f - tmp_cpc2.tsv \\
        | awk 'BEGIN{FS=OFS="\\t"}{print \$1, "0.001"}' \\
        >> tmp_tpm.tsv

    # Re-sort the augmented TPM file
    sort tmp_tpm.tsv -o tmp_tpm_sorted.tsv

    # Sort CPC2 file to align with TPM
    sort tmp_cpc2.tsv -o tmp_cpc2_sorted.tsv

    # Join on transcript ID and emit header + data
    paste tmp_cpc2_sorted.tsv tmp_tpm_sorted.tsv \\
        | bedtools groupby -g 1 -c 2,4 -o collapse \\
        | sed 's/0\\.000000/0.000001/g' \\
        | awk 'BEGIN{print "tid\\tCPC\\ttpm"; FS=OFS="\\t"}{print}' \\
        > tmp_cpc2_tpm.tsv

    # Merge BUSCO scores column if available
    if ${has_busco}; then
        # Build a busco_score lookup (skip header), default 0 for missing transcripts
        tail -n +2 ${busco_scores_tsv} | sort > tmp_busco_sorted.tsv
        # Left-join: for each transcript in CPC2/TPM, add busco_score (0 if not in busco file)
        awk 'BEGIN{FS=OFS="\\t"}
            NR==FNR { busco[\$1]=\$2; next }
            FNR==1  { print \$0, "busco_score"; next }
            { print \$0, ((\$1 in busco) ? busco[\$1] : "0") }
        ' tmp_busco_sorted.tsv tmp_cpc2_tpm.tsv > external_scores_tpm_cpc2.tsv
    else
        mv tmp_cpc2_tpm.tsv external_scores_tpm_cpc2.tsv
    fi
    """

    stub:
    """
    touch external_scores_tpm_cpc2.tsv
    """
}
