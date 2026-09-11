process RUN_DBCAN {
    tag "run_dbcan"
    label 'process_high'
    container 'ghcr.io/bcb-unl/run_dbcan:5.2.9'

    publishDir "${params.outdir}/functional_annotation/dbcan", mode: 'copy',
        saveAs: { it.replaceFirst(/^dbcan_out\//, '') },
        pattern: "dbcan_out/*"

    input:
    path proteins
    path gff
    path db_dir, stageAs: 'dbcan_db'

    output:
    path "dbcan_out/overview.tsv",         emit: overview
    path "dbcan_out/cgc_standard_out.tsv", emit: cgc_standard
    path "dbcan_out/cgc.gff",              emit: cgc_gff
    path "dbcan_out/*"

    script:
    """
    # Add protein_id to CDS features so NCBI_euk parser can map genes → proteins
    python3 - <<'PYGFF'
import re
with open('${gff}') as fin, open('annot_protein_id.gff', 'w') as fout:
    for line in fin:
        if line.startswith('#') or not line.strip():
            fout.write(line)
            continue
        parts = line.rstrip('\\n').split('\\t')
        if len(parts) >= 9 and parts[2] == 'CDS':
            m = re.search(r'Parent=([^;]+)', parts[8])
            if m:
                parts[8] = parts[8].rstrip(';') + ';protein_id=' + m.group(1)
        fout.write('\\t'.join(parts) + '\\n')
PYGFF

    run_dbcan easy_CGC \\
        --mode protein \\
        --input_raw_data ${proteins} \\
        --db_dir dbcan_db \\
        --output_dir dbcan_out \\
        --gff_type NCBI_euk \\
        --input_gff annot_protein_id.gff \\
        --threads ${task.cpus} \\
        --methods diamond \\
        --methods hmm \\
        --methods dbCANsub
    """

    stub:
    """
    mkdir -p dbcan_out
    printf 'Gene ID\\tEC#\\tdbCAN_hmm\\tdbCAN_sub\\tDIAMOND\\t#ofTools\\tRecommend Results\\n' > dbcan_out/overview.tsv
    printf 'CGC#\\tGene Type\\tContig ID\\tProtein ID\\tGene Start\\tGene Stop\\tGene Strand\\tGene Annotation\\n' > dbcan_out/cgc_standard_out.tsv
    touch dbcan_out/cgc.gff
    """
}
