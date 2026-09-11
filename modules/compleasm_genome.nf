process COMPLEASM_GENOME {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/compleasm:0.2.8--pyh106432d_0'

    publishDir "${params.outdir}/compleasm_genome", mode: 'copy', pattern: "compleasm_out/summary.txt"
    publishDir "${params.outdir}/compleasm_genome", mode: 'copy', pattern: "busco.gff"

    input:
    tuple val(meta), path(fasta)
    val   lineage

    output:
    tuple val(meta), path("compleasm_out/summary.txt"), emit: summary
    tuple val(meta), path("busco.gff"),                 emit: busco_anno_gff
    path "mb_downloads",                                emit: busco_db

    script:
    """
    compleasm run \\
        -a ${fasta} \\
        -o compleasm_out \\
        -l ${lineage}_odb12 \\
        -t ${task.cpus}

    # compleasm resolves odb10 lineages to odb12 equivalents; rename to expected path
    expected="compleasm_out/${lineage}_odb12"
    if [ ! -d "\$expected" ]; then
        actual=\$(find compleasm_out -maxdepth 1 -mindepth 1 -type d | head -1)
        [ -n "\$actual" ] && mv "\$actual" "\$expected"
    fi

    awk '!/^#/' "\$expected/miniprot_output.gff" | sed 's/ /_/g' > busco.gff
    """

    stub:
    """
    mkdir -p compleasm_out/${lineage}_odb12
    mkdir -p mb_downloads/${lineage}_odb12
    touch compleasm_out/summary.txt
    touch busco.gff
    """
}
