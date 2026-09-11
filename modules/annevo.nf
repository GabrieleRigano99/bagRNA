process ANNEVO {
    tag "annevo"
    label 'process_high'
    container 'gabrielerigano/annevo:v2.2.2'
    containerOptions "-v \$PWD:/data --shm-size=${params.ram_annevo}${params.use_gpu ? ' --gpus all' : ''}"

    publishDir "${params.outdir}/annevo", mode: 'copy'

    input:
    path  fasta
    val   species
    val   strain
    val   lineage

    output:
    path "ANNEVO.gff", emit: annevo_gff

    script:
    def model_map = [
        'fungi'       : 'ANNEVO_Fungi.pt',
        'land_plant'  : 'ANNEVO_Embryophyta.pt',
        'vertebrate'  : 'ANNEVO_Vertebrate_other.pt',
        'invertebrate': 'ANNEVO_Invertebrate.pt',
        'magnoliopsida': 'ANNEVO_Magnoliopsida.pt',
        'mammalia'    : 'ANNEVO_Mammalia2.pt',
        'insecta'     : 'ANNEVO_Insecta.pt'
    ]
    def model = model_map[lineage]
    if (!model) error "Unknown --annevo_lineage '${lineage}'. Valid values: ${model_map.keySet().join(', ')}"
    """
    python /opt/ANNEVO/annotation.py \\
        --genome /data/${fasta} \\
        --model_path /opt/ANNEVO/ANNEVO_model/${model} \\
        --output /data/ANNEVO.gff
    """

    stub:
    """
    touch ANNEVO.gff
    """
}
