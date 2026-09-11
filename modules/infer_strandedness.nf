process INFER_STRANDEDNESS {
    tag "infer_strandedness"
    label 'process_low'
    container 'quay.io/biocontainers/rseqc:5.0.4--pyhdfd78af_1'

    publishDir "${params.outdir}/strandedness", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf

    output:
    path "infer_experiment.txt",  emit: infer_txt
    path "strandedness.txt",      emit: strandedness_file
    path "orientation.txt",       emit: orientation_file
    path "trinity_lib_type.txt",  emit: trinity_lib_type_file

    script:
    """
    python3 << 'PYEOF'
    from collections import defaultdict

    # ── GTF → BED12 ──────────────────────────────────────────────────────────
    transcripts = {}

    with open("${gtf}") as fh:
        for line in fh:
            if line.startswith('#'): continue
            fields = line.rstrip('\\n').split('\\t')
            if len(fields) < 9: continue
            chrom, _, feature, start, end, _, strand, _, attrs = fields
            if feature not in ('transcript', 'exon'): continue

            tid = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    tid = attr.split('"')[1] if '"' in attr else attr.split(' ')[-1]
                    break
            if not tid: continue

            start, end = int(start) - 1, int(end)

            if feature == 'transcript':
                transcripts[tid] = {'chrom': chrom, 'strand': strand,
                                    'start': start, 'end': end, 'exons': []}
            elif feature == 'exon':
                if tid not in transcripts:
                    transcripts[tid] = {'chrom': chrom, 'strand': strand,
                                        'start': start, 'end': end, 'exons': []}
                transcripts[tid]['exons'].append((start, end))

    with open('annevo.bed', 'w') as out:
        for tid, info in transcripts.items():
            exons = sorted(info['exons'])
            if not exons: continue
            t_start = min(e[0] for e in exons)
            t_end   = max(e[1] for e in exons)
            sizes   = ','.join(str(e[1] - e[0]) for e in exons) + ','
            starts  = ','.join(str(e[0] - t_start) for e in exons) + ','
            out.write(f"{info['chrom']}\\t{t_start}\\t{t_end}\\t{tid}\\t0\\t"
                      f"{info['strand']}\\t{t_start}\\t{t_end}\\t0\\t"
                      f"{len(exons)}\\t{sizes}\\t{starts}\\n")
    PYEOF

    infer_experiment.py -r annevo.bed -i ${bam} -s 200000 > infer_experiment.txt

    python3 << 'PYEOF'
    import re

    with open('infer_experiment.txt') as fh:
        text = fh.read()

    # Paired-end: parse both fractions
    m1 = re.search(r'Fraction.*?"1\\+\\+,1--,2\\+-,2-\\+".*?:\\s*([\\d.]+)', text)
    m2 = re.search(r'Fraction.*?"1\\+-,1-\\+,2\\+\\+,2--".*?:\\s*([\\d.]+)', text)

    threshold = 0.6

    if m1 and m2:
        f_sense    = float(m1.group(1))   # 1++,1--,2+-,2-+ → read1 sense → secondstrand
        f_antisense = float(m2.group(1))  # 1+-,1-+,2++,2-- → read2 sense → firststrand (dUTP)
        if f_sense >= threshold:
            strandedness    = 'secondstrand'
            orientation     = 'FR'
            trinity_lib     = 'FR'
        elif f_antisense >= threshold:
            strandedness    = 'firststrand'
            orientation     = 'FR'
            trinity_lib     = 'RF'
        else:
            strandedness    = 'unstranded'
            orientation     = 'unstranded'
            trinity_lib     = ''
    else:
        # Single-end or unparseable → fall back to defaults
        strandedness = 'secondstrand'
        orientation  = 'FR'
        trinity_lib  = 'FR'

    with open('strandedness.txt',     'w') as f: f.write(strandedness)
    with open('orientation.txt',      'w') as f: f.write(orientation)
    with open('trinity_lib_type.txt', 'w') as f: f.write(trinity_lib)

    print(f"[INFER_STRANDEDNESS] strandedness={strandedness}  orientation={orientation}  trinity_lib={trinity_lib}")
    print(f"  fractions → sense: {m1.group(1) if m1 else 'n/a'}  antisense: {m2.group(1) if m2 else 'n/a'}")
    PYEOF
    """

    stub:
    """
    echo 'secondstrand' > strandedness.txt
    echo 'FR'           > orientation.txt
    echo 'FR'           > trinity_lib_type.txt
    touch infer_experiment.txt
    """
}
