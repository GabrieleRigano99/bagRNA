process SIGNALP6 {
    tag "signalp6"
    label 'process_high'
    container 'interpro/signalp:6.0i'
    containerOptions {
        def pkg = file(params.signalp_path).toAbsolutePath()
        "${params.use_gpu ? '--gpus all' : ''} -v ${pkg}:/tools".trim()
    }

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path proteins_faa

    output:
    path "signalp_output/prediction_results.txt", emit: signalp_results

    script:
    def convert_cmd = params.use_gpu ? 'signalp6_convert_models gpu' : ''
    """
    export HOME=/tmp
    pip install /tools/signalp-6-package --user --quiet --no-cache-dir
    export PATH="\${HOME}/.local/bin:\${PATH}"
    ${convert_cmd}
    signalp6 \\
        --fastafile ${proteins_faa} \\
        --organism eukarya \\
        --format txt \\
        --output_dir signalp_output \\
        --torch_num_threads ${task.cpus} \\
        --model_dir /tools/signalp-6-package/models/
    """

    stub:
    """
    mkdir -p signalp_output
    touch signalp_output/prediction_results.txt
    """
}
