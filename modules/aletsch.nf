process ALETSCH {
    tag "aletsch_assembly"
    label 'process_high'
    container 'gabrielerigano/aletsch:1.1.3'

    publishDir "${params.outdir}/transcript_assembly", mode: 'copy'

    input:
    path bams

    output:
    path "aletsch.gtf", emit: aletsch_gtf

    script:
    def bam_list = bams instanceof List ? bams : [bams]
    """
    # Write the input BAM list required by aletsch
    for bam in ${bam_list.join(' ')}; do
        samtools index -@ ${task.cpus} "\$bam" 2>/dev/null || true
        echo "\$bam \${bam}.bai paired_end"
    done > input_bam_list.txt

    mkdir -p profile_dir gtf_dir

    # Profile pass
    aletsch --profile \\
        -i input_bam_list.txt \\
        -p profile_dir \\
        > preprocess_aletsch.log 2>&1

    # Assembly pass
    aletsch \\
        -i input_bam_list.txt \\
        -o raw_aletsch.gtf \\
        -p profile_dir \\
        --output_gtf_dir gtf_dir \\
        --max_threads ${task.cpus} \\
        --boost_precision \\
        --min_splice_bundary_hits 5 \\
        --min_transcript_length_base 200 \\
        --min_mapping_quality 20 \\
        --output_single_exon_transcripts \\
        > aletsch_assembly.log 2>&1

    # Merge with stringtie to produce a clean GTF
    stringtie --merge -o aletsch.gtf raw_aletsch.gtf
    """

    stub:
    """
    touch aletsch.gtf
    """
}
