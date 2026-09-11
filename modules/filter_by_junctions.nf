process FILTER_BY_JUNCTIONS {
    tag "filter_by_junctions"
    label 'process_low'
    container 'python:3.11-slim'

    publishDir "${params.outdir}/transcript_assembly", mode: 'copy',
        saveAs: { fn -> fn.startsWith('validated_') ? fn : null }

    input:
    path(stringtie_gtf, stageAs: 'stringtie.gtf')
    path(aletsch_gtf,   stageAs: 'aletsch.gtf')
    path(trinity_gff,   stageAs: 'trinity.gff')
    path portcullis_bed

    output:
    path "validated_stringtie.gtf", emit: stringtie_gtf
    path "validated_aletsch.gtf",   emit: aletsch_gtf
    path "validated_trinity.gff",   emit: trinity_gff

    script:
    """
    if [ -s trinity.gff ]; then
        filter_by_junctions.py \\
            --junctions ${portcullis_bed} \\
            --input  stringtie.gtf aletsch.gtf trinity.gff \\
            --output validated_stringtie.gtf validated_aletsch.gtf validated_trinity.gff
    else
        filter_by_junctions.py \\
            --junctions ${portcullis_bed} \\
            --input  stringtie.gtf aletsch.gtf \\
            --output validated_stringtie.gtf validated_aletsch.gtf
        touch validated_trinity.gff
    fi
    """

    stub:
    """
    touch validated_stringtie.gtf validated_aletsch.gtf validated_trinity.gff
    """
}
