process MIKADO_PREPARE {
    tag "mikado_prepare"
    label 'process_medium'
    container 'baderlab/mikado:ubuntu22_mikado2.3.2'

    publishDir "${params.outdir}/mikado",              mode: 'copy'
    publishDir "${params.outdir}/mikado", mode: 'copy', pattern: "mikado_op/mikado_prepared.gtf",   saveAs: { "mikado_prepared.gtf" }
    publishDir "${params.outdir}/mikado", mode: 'copy', pattern: "mikado_op/mikado_prepared.fasta", saveAs: { "mikado_prepared.fasta" }

    input:
    path fasta
    path junctions_bed
    path prot_evidence
    val  scoring_yaml
    path(scoring_file,            stageAs: 'scoring_custom.yaml')
    path(helixer_gff,             stageAs: 'helixer.gff')
    path(annevo_gtf,              stageAs: 'annevo.gtf')
    path(stringtie_gtf,           stageAs: 'stringtie.gtf')
    path(trinity_gff,             stageAs: 'trinity.gff')
    path(aletsch_gtf,             stageAs: 'aletsch.gtf')
    path(busco_gff,               stageAs: 'busco.gff')
    path(miniprot_gtf,            stageAs: 'miniprot.gtf')
    path(liftover_gff,            stageAs: 'liftover.gff')
    path(transcript_evidence_gff, stageAs: 'transcript_evidence.gff')

    output:
    path "mikado_op/",                      emit: mikado_dir
    path "mikado_op/mikado_prepared.fasta", emit: prepared_fasta
    path "mikado_op/mikado_prepared.gtf",   emit: prepared_gtf
    path "configuration.yaml",              emit: configuration

    script:
    """
    mkdir -p mikado_op

    if [ -s "scoring_custom.yaml" ]; then
        cp scoring_custom.yaml ${scoring_yaml}
    else
        cp /usr/local/lib/python3.10/dist-packages/Mikado/configuration/scoring_files/HISTORIC/${scoring_yaml} ./
    fi

    # Build the Mikado list dynamically — only include sources that are present
    # Columns: filename  label  strand_specific  score  is_reference  exclude_redundant  strip_cds
    > mikado_list.tsv

    # Always-present sources (pipeline always produces these)
    printf "stringtie.gtf\\tStringtie\\tTrue\\t1\\tFalse\\tTrue\\tTrue\\n"       >> mikado_list.tsv
    printf "aletsch.gtf\\tAletsch\\tTrue\\t1\\tFalse\\tTrue\\tTrue\\n"           >> mikado_list.tsv
    printf "busco.gff\\tBusco\\tTrue\\t2\\tFalse\\tTrue\\tTrue\\n"               >> mikado_list.tsv
    printf "miniprot.gtf\\tProteinEvidence\\tTrue\\t3\\tFalse\\tTrue\\tTrue\\n"  >> mikado_list.tsv

    # Optional sources — only added when the staged file is non-empty
    [ -s helixer.gff ]              && printf "helixer.gff\\tHelixer\\tTrue\\t3\\tTrue\\tTrue\\tTrue\\n"                         >> mikado_list.tsv
    [ -s annevo.gtf ]               && printf "annevo.gtf\\tAnnevo\\tTrue\\t4\\tTrue\\tTrue\\tTrue\\n"                           >> mikado_list.tsv
    [ -s trinity.gff ]              && printf "trinity.gff\\tTrinity\\tTrue\\t-0.5\\tFalse\\tTrue\\tTrue\\n"                     >> mikado_list.tsv
    [ -s liftover.gff ]             && printf "liftover.gff\\tLiftover\\tTrue\\t2\\tFalse\\tTrue\\tTrue\\n"                      >> mikado_list.tsv
    [ -s transcript_evidence.gff ]  && printf "transcript_evidence.gff\\tTranscriptEvidence\\tTrue\\t3\\tFalse\\tTrue\\tTrue\\n" >> mikado_list.tsv

    mikado configure \\
        --list mikado_list.tsv \\
        --reference ${fasta} \\
        --junctions ${junctions_bed} \\
        --max-intron-length ${params.max_intron_length} \\
        --strand-specific \\
        --mode permissive \\
        --check-references \\
        --use-transdecoder \\
        --scoring ${scoring_yaml} \\
        --blast_targets ${prot_evidence} \\
        --threads ${task.cpus} \\
        --out-dir mikado_op \\
        --yaml configuration.yaml

    mikado prepare \\
        --proc ${task.cpus} \\
        --seed 0 \\
        --json-conf configuration.yaml
    """

    stub:
    """
    mkdir -p mikado_op
    touch mikado_op/mikado_prepared.fasta
    touch mikado_op/mikado_prepared.gtf
    touch configuration.yaml
    """
}
