// Miniprot aligns every protein in the evidence set independently, so a
// multi-species reference proteome (e.g. a full BUSCO/compleasm lineage
// protein set) produces massive paralog/ortholog redundancy at each locus
// — thousands of near-identical hits from different reference species all
// landing on the same gene. AGAT_MERGE_ABINITIO's overlap-based union has
// no way to tell those apart from real alternative isoforms, so it inflates
// single loci into dozens/hundreds of "transcripts". gffcompare's structural
// comparison (identical/contained intron-chain clustering) collapses that
// redundancy properly; we use it here only to pick one representative
// Miniprot transcript per unique structure, then hand the reduced set to
// AGAT_MERGE_ABINITIO, which still does the real gene/CDS model building.
process GFFCOMPARE_DEDUP_MINIPROT {
    tag "gffcompare_dedup_miniprot"
    label 'process_medium'
    container 'quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_1'

    input:
    path helixer_gff
    path annevo_gtf
    path miniprot_gtf
    path barrnap_gff

    output:
    path "miniprot_dedup.gtf", emit: miniprot_dedup_gtf

    script:
    """
    ARGS=""
    for f in ${helixer_gff} ${annevo_gtf} ${miniprot_gtf} ${barrnap_gff}; do
        [ -s "\$f" ] && ARGS="\$ARGS \$f"
    done

    gffcompare -o gffcmp \$ARGS

    if [ -s gffcmp.tracking ]; then
        # Column 7 is the miniprot (3rd query file) slot; take its oId
        # (2nd '|'-delimited field) for every locus where miniprot
        # contributed the retained representative structure.
        awk -F'\\t' '\$7!="-"{split(\$7,a,"|"); print a[2]}' gffcmp.tracking \\
            | sort -u \\
            | sed 's/^/transcript_id "/; s/\$/"/' \\
            > keep_patterns.txt

        if [ -s keep_patterns.txt ]; then
            grep -F -f keep_patterns.txt ${miniprot_gtf} > miniprot_dedup.gtf || true
        else
            : > miniprot_dedup.gtf
        fi
    else
        : > miniprot_dedup.gtf
    fi
    """

    stub:
    """
    touch miniprot_dedup.gtf
    """
}
