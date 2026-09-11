process STAR_ALIGN {
    tag "star_align"
    label 'process_high'
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_1'

    publishDir "${params.outdir}/star_align", mode: 'copy', pattern: "*.log"
    publishDir "${params.outdir}/star_align", mode: 'copy', pattern: "*.bam"

    input:
    path  manifest
    path  star_index
    path  reads  // staged to trigger Docker volume-mount of the reads directory

    output:
    path "star_Aligned.sortedByCoord.out.bam", emit: bam
    path "star_Log.final.out",                 emit: log_final

    script:
    // Use 90 % of the task's allocated memory for BAM sorting; user can override via --limitBAMsortRAM
    def bam_sort_ram = params.limitBAMsortRAM ?: (long)(task.memory.toBytes() * 0.9)
    """
    # Build 3-column manifest required by STAR: R1<TAB>R2<TAB>ID:sample<TAB>SM:sample
    # Attributes within the RG column must be tab-separated so STAR parses ID correctly
    awk 'BEGIN{OFS="\\t"}{
        r2 = (NF>=2 && \$2!="" && \$2!="-") ? \$2 : "-"
        id = (NF>=3 && \$3!="") ? \$3 : "sample"NR
        rg = "ID:"id"\\tSM:"id
        print \$1, r2, rg
    }' ${manifest} > star_rg_manifest.tsv

    # --outSAMattrRGline: attributes tab-separated within each group, groups comma-separated
    rg_line=\$(awk '{ id=(NF>=3 && \$3!="") ? \$3 : "sample"NR
                     printf "%sID:%s\\tSM:%s", (NR>1 ? " , " : ""), id, id }' ${manifest})

    STAR \\
        --runThreadN ${task.cpus} \\
        --twopassMode Basic \\
        --genomeDir ${star_index} \\
        --genomeSAindexNbases ${params.genomeSAindexNbases} \\
        --readFilesManifest star_rg_manifest.tsv \\
        --readFilesCommand zcat \\
        --outSAMattrRGline \${rg_line} \\
        --outFileNamePrefix star_ \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS NM MD RG \\
        --outSAMstrandField intronMotif \\
        --alignIntronMax ${params.max_intron_length} \\
        --limitBAMsortRAM ${bam_sort_ram} \\
        --outBAMcompression 10
    """

    stub:
    """
    touch star_Aligned.sortedByCoord.out.bam
    touch star_Log.final.out
    """
}
