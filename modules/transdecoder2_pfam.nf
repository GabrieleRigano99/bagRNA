process TRANSDECODER2_PFAM_SCAN {
    tag "td2_pfam"
    label 'process_high'
    container 'ghcr.io/bcb-unl/run_dbcan:5.2.9'

    input:
    path pep
    path(pfam_hmm, stageAs: 'pfam_src.hmm')

    output:
    path "pfam.domtblout", emit: domtblout

    script:
    """
    python3 - << 'PYEOF'
import pyhmmer

with pyhmmer.easel.SequenceFile("${pep}", digital=True) as sf:
    sequences = sf.read_block()

with pyhmmer.plan7.HMMFile("pfam_src.hmm") as hf:
    with open("pfam.domtblout", "wb") as out:
        for hits in pyhmmer.hmmer.hmmscan(sequences, hf, cpus=${task.cpus}, E=1e-5):
            hits.write(out, format="domains")
PYEOF
    """

    stub:
    """
    touch pfam.domtblout
    """
}
