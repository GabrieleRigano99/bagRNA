// table2asn.nf
// Structural annotation → NCBI submission files.
// Used at the end of structural annotation (step 24) to produce .gbf files
// that feed into AntiSMASH.  For fully-annotated submission files (with
// product names, EC numbers, GO dbxrefs) see modules/ncbi_submission.nf.

process TABLE2ASN {
    tag "table2asn"
    label 'process_medium'
    container 'staphb/ncbi-table2asn'

    publishDir "${params.outdir}/table2asn", mode: 'copy'

    input:
    path final_gff
    path fasta
    path submission_template
    val  species
    val  strain
    val  codon_table
    val  locus_tag

    output:
    path "table2asn/",       emit: table2asn_dir
    path "table2asn/*.gbf",  emit: gbf_files

    script:
    """
    mkdir -p table2asn

    # No -Z (discrepancy report): reproducibly hangs indefinitely at this stage
    # on this pipeline's raw (pre-functional-annotation) structural GFF — twice,
    # both times stalling right after .gbf/.sqn were written, never producing
    # .dr/.stats (2026-08-28). Nothing downstream consumes .dr/.stats (only
    # .gbf, by AntiSMASH), so dropping it removes the hang with no pipeline
    # impact. NCBI_SUBMISSION's separate table2asn call (on the fully
    # functionally-annotated GFF) keeps -Z — it has run successfully every time.
    table2asn \\
        -V vb \\
        -M n \\
        -J \\
        -c ewf \\
        -euk \\
        -t ${submission_template} \\
        -gaps-min 10 \\
        -l paired-ends \\
        -locus-tag-prefix ${locus_tag} \\
        -j "[organism=${species}] [strain=${strain}] [gcode=${codon_table}]" \\
        -i ${fasta} \\
        -f ${final_gff} \\
        -outdir table2asn \\
        -verbose
    """

    stub:
    """
    mkdir -p table2asn
    touch table2asn/stub_output.gbf
    """
}
