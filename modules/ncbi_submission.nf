// ncbi_submission.nf
// Generate fully-annotated NCBI submission files (.sqn + .gbf) for genome deposition.
// Runs after ANNOTATE_FUNCTIONAL: uses gff3_to_tbl.py to merge GFF3 structure with
// bagRNA's functional_annotation.tsv (product names, EC numbers, GO db_xrefs, notes)
// into a NCBI feature table (.tbl), then calls table2asn to produce .sqn and .gbf.

process GFF3_TO_TBL {
    tag "gff3_to_tbl"
    label 'process_low'
    container 'python:3.11-slim'
    containerOptions '-v /usr/bin/ps:/usr/bin/ps'

    input:
    path final_gff
    path annotation_tsv
    val  species
    val  strain

    output:
    path "${prefix}.tbl", emit: tbl

    script:
    prefix = "${species}_${strain}".replaceAll(/\s+/, '_')
    """
    python3 ${workflow.projectDir}/bin/gff3_to_tbl.py \\
        ${final_gff} \\
        ${annotation_tsv} \\
        ${prefix}.tbl
    """

    stub:
    prefix = "${species}_${strain}".replaceAll(/\s+/, '_')
    """
    touch ${prefix}.tbl
    """
}

process NCBI_SUBMISSION {
    tag "ncbi_submission"
    label 'process_medium'
    container 'staphb/ncbi-table2asn'

    publishDir "${params.outdir}/ncbi_submission", mode: 'copy'

    input:
    path fasta               // genome FASTA
    path tbl                 // feature table from GFF3_TO_TBL
    path submission_template // NCBI .sbt template (or NO_FILE)
    val  species
    val  strain
    val  codon_table
    val  locus_tag

    output:
    path "${prefix}.sqn",   emit: sqn,   optional: true
    path "${prefix}.gbf",   emit: gbf,   optional: true
    path "${prefix}.tbl",   emit: tbl
    path "${prefix}.val",   emit: val,   optional: true
    path "${prefix}.stats", emit: stats, optional: true

    script:
    prefix  = "${species}_${strain}".replaceAll(/\s+/, '_')
    sbt_arg = (submission_template.name != 'NO_FILE') ? "-t ${submission_template}" : ''
    """
    # Generate .sqn and .gbf submission files
    table2asn \\
        -V vb \\
        -M n \\
        -J \\
        -c fx \\
        -euk \\
        ${sbt_arg} \\
        -gaps-min 10 \\
        -l paired-ends \\
        -locus-tag-prefix ${locus_tag} \\
        -j "[organism=${species}] [strain=${strain}] [gcode=${codon_table}]" \\
        -i ${fasta} \\
        -f ${prefix}.tbl \\
        -o ${prefix}.sqn \\
        -w ${prefix}.stats \\
        -Z

    # Normalise output extension (.gbk → .gbf differs across table2asn versions)
    if [ -f ${prefix}.gbk ] && [ ! -f ${prefix}.gbf ]; then
        mv ${prefix}.gbk ${prefix}.gbf
    fi
    """

    stub:
    prefix = "${species}_${strain}".replaceAll(/\s+/, '_')
    """
    touch ${prefix}.sqn ${prefix}.gbf ${prefix}.tbl ${prefix}.val ${prefix}.stats
    """
}
