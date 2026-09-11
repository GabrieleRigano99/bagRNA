process CORRECT_MISORIENTED {
    tag "correct_misoriented"
    label 'process_low'
    container 'python:3.11-slim'

    publishDir "${params.outdir}/transcript_assembly", mode: 'copy',
        saveAs: { fn -> fn.startsWith('corrected_') ? fn : null }

    input:
    path(stringtie_gtf, stageAs: 'stringtie.gtf')
    path(aletsch_gtf,   stageAs: 'aletsch.gtf')
    path(trinity_gff,   stageAs: 'trinity.gff')
    path miniprot_gtf

    output:
    path "corrected_stringtie.gtf", emit: stringtie_gtf
    path "corrected_aletsch.gtf",   emit: aletsch_gtf
    path "corrected_trinity.gff",   emit: trinity_gff

    script:
    """
    if [ -s trinity.gff ]; then
        correct_misoriented.py \\
            --miniprot ${miniprot_gtf} \\
            --input  stringtie.gtf aletsch.gtf trinity.gff \\
            --output corrected_stringtie.gtf corrected_aletsch.gtf corrected_trinity.gff
    else
        correct_misoriented.py \\
            --miniprot ${miniprot_gtf} \\
            --input  stringtie.gtf aletsch.gtf \\
            --output corrected_stringtie.gtf corrected_aletsch.gtf
        touch corrected_trinity.gff
    fi
    """

    stub:
    """
    touch corrected_stringtie.gtf corrected_aletsch.gtf corrected_trinity.gff
    """
}
