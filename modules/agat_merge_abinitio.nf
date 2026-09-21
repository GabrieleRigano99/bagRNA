// No-RNAseq structural annotation: merges Helixer + ANNEVO (ab initio) +
// Miniprot (protein-to-genome) + Barrnap + tRNAscan-SE into one GFF3 with
// AGAT's own parser/merger, rather than Mikado's evidence-weighted
// consensus (which needs RNA-seq-derived junctions/TPM to score against).
// Any input that's the NO_FILE sentinel (empty, e.g. --no_helixer) or a
// genuinely empty prediction (no rRNA/tRNA found) is skipped.
process AGAT_MERGE_ABINITIO {
    tag "agat_merge_abinitio"
    label 'process_medium'
    container 'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1'

    input:
    path helixer_gff
    path annevo_gff
    path miniprot_gtf
    path barrnap_gff
    path trnascan_gff

    output:
    path "merged_abinitio.gff", emit: merged_gff

    script:
    """
    ARGS=""
    for f in ${helixer_gff} ${annevo_gff} ${miniprot_gtf} ${barrnap_gff} ${trnascan_gff}; do
        [ -s "\$f" ] && ARGS="\$ARGS --gff \$f"
    done

    agat_sp_merge_annotations.pl \$ARGS -o merged_abinitio.gff

    # AGAT doesn't recognize tRNAscan-SE's "pseudogene" gene_biotype (used
    # for low-confidence tRNA hits) and falls back to a generic "RNA"
    # transcript type for its child feature. table2asn rejects "RNA" as an
    # invalid feature type outright. These are always tRNA pseudogenes
    # (nothing else in this pipeline's inputs produces bare "RNA"), so
    # relabel them tRNA.
    awk -F'\\t' 'BEGIN{OFS="\\t"} \$3=="RNA"{\$3="tRNA"} {print}' merged_abinitio.gff > merged_abinitio.fixed.gff
    mv merged_abinitio.fixed.gff merged_abinitio.gff
    """

    stub:
    """
    touch merged_abinitio.gff
    """
}
