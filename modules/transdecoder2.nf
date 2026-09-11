process TRANSDECODER2_LONGORFS {
    tag "td2_longorfs"
    label 'process_high'
    label 'use_gpu'
    container 'gabrielerigano/td2:1.1.0-cuda12.4'

    input:
    path prepared_fasta
    val  strandedness

    output:
    path "transdecoder",                    emit: td2_dir
    path "transdecoder/longest_orfs.pep",   emit: pep

    script:
    def strand_flag = (strandedness && strandedness != 'unstranded') ? '-S' : ''
    """
    TD2.LongOrfs \\
        -G ${params.codon_table} \\
        --threads ${task.cpus} \\
        --complete-orfs-only \\
        ${strand_flag} \\
        -O transdecoder \\
        -t ${prepared_fasta}
    """

    stub:
    """
    mkdir -p transdecoder
    touch transdecoder/longest_orfs.pep
    """
}

process TRANSDECODER2_PREDICT {
    tag "td2_predict"
    label 'process_high'
    label 'use_gpu'
    container 'gabrielerigano/td2:1.1.0-cuda12.4'
    cache 'deep'

    publishDir "${params.outdir}/transdecoder", mode: 'copy'

    input:
    path prepared_fasta
    path(td2_dir,     stageAs: 'transdecoder')
    path hmmer_hits   // pfam domtblout; NO_FILE to skip
    path blastp_hits  // outfmt6 hits; NO_FILE to skip

    output:
    path "${prepared_fasta.name}.TD2.bed",  emit: bed
    path "${prepared_fasta.name}.TD2.gff3", emit: gff3

    script:
    def hmmer_arg  = hmmer_hits.name  != 'NO_FILE' ? "--retain-hmmer_hits ${hmmer_hits}"  : ''
    def blastp_arg = blastp_hits.name != 'NO_FILE' ? "--retain-blastp_hits ${blastp_hits}" : ''
    """
    TD2.Predict \\
        -G ${params.codon_table} \\
        -O transdecoder \\
        -t ${prepared_fasta} \\
        ${hmmer_arg} \\
        ${blastp_arg}

    # Fix minus-strand ORFs with empty thickEnd (TD2 bug)
    python3 - << 'PYEOF'
import re, os

bed = "${prepared_fasta.name}.TD2.bed"
fixed = bed + ".fix"
with open(bed) as fin, open(fixed, "w") as fout:
    for line in fin:
        if line.startswith("track") or "\\t" not in line:
            fout.write(line)
            continue
        f = line.rstrip("\\n").split("\\t")
        if len(f) >= 8 and f[7] == "":
            m = re.search(r':(\\d+)-(\\d+)\\([+-]\\)\$', f[3])
            if m:
                f[7] = m.group(2)
        fout.write("\\t".join(f) + "\\n")

os.replace(fixed, bed)
PYEOF
    """

    stub:
    """
    touch "${prepared_fasta.name}.TD2.bed"
    touch "${prepared_fasta.name}.TD2.gff3"
    """
}
