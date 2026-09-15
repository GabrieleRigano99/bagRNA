// download_databases.nf
// bagRNA database setup workflow
// Run with: nextflow run main.nf -entry SETUP --db_dir /path/to/databases
// All databases are placed under db_dir/:
//
//   db_dir/
//   ├── eggnog_data/                             (→ --eggnog_data_dir db_dir/eggnog_data)
//   ├── Rfam.cm, Rfam.clanin                     (→ --databases db_dir)
//   ├── dbcan/                                   (→ --dbcan_db db_dir/dbcan)
//   ├── merops/merops_pepunit.dmnd               (→ --merops_db db_dir/merops/merops_pepunit.dmnd)
//   └── phi_base/phi_base.fas                    (→ --phi_base_db db_dir/phi_base/phi_base.fas)

nextflow.enable.dsl = 2

// ── 1. eggNOG-mapper databases (v7) ──────────────────────────────────────────
// ~43 GB; eggnog.db, eggnog_proteins.dmnd, go-basic.obo, taxonomy files.
// Mirrored from the eggnog-mapper 3.0 (v7 DB) data release — the bundled
// download_eggnog_data.py only ever fetches v5.0.2, so this pulls the v7
// tree directly instead. Publishes as db_dir/eggnog_data/; pass that whole
// directory to --eggnog_data_dir (matches what modules/eggnog.nf, running
// gabrielerigano/eggnog-mapper:3.0.0-beta6, expects).
process DOWNLOAD_EGGNOG {
    tag "eggnog"
    label 'process_long'
    container 'quay.io/biocontainers/wget:1.20.1'

    publishDir "${params.db_dir}", mode: 'copy'

    output:
    path "eggnog_data"

    script:
    """
    mkdir -p eggnog_data
    wget -q -r -np -nH --cut-dirs=3 -X mmseqs -R "index.html*" -P eggnog_data \\
        https://data.cgmlab.org/eggnog-mapper/emapper-3.0/data/
    """

    stub:
    """
    mkdir -p eggnog_data
    touch eggnog_data/eggnog.db eggnog_data/eggnog_proteins.dmnd
    """
}

// ── 2. Rfam covariance models (for Infernal ncRNA search) ────────────────────
// ~1 GB Rfam.cm; cmpress indexes it for fast cmsearch
process DOWNLOAD_RFAM {
    tag "rfam"
    label 'process_medium'
    container 'quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_3'

    publishDir "${params.db_dir}", mode: 'copy'

    output:
    path "Rfam.cm"
    path "Rfam.cm.i1f"
    path "Rfam.cm.i1i"
    path "Rfam.cm.i1m"
    path "Rfam.cm.i1p"
    path "Rfam.clanin"

    script:
    """
    wget -q https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
    wget -q https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
    gunzip Rfam.cm.gz
    cmpress Rfam.cm
    """

    stub:
    """
    touch Rfam.cm Rfam.cm.i1f Rfam.cm.i1i Rfam.cm.i1m Rfam.cm.i1p Rfam.clanin
    """
}

// ── 3. dbCAN databases (CAZyme + CGC annotation) ─────────────────────────────
// CAZy.dmnd, dbCAN.hmm, dbCAN-sub.hmm, tcdb.dmnd, tf HMMs, stp.dmnd
process DOWNLOAD_DBCAN {
    tag "dbcan"
    label 'process_high'
    container 'quay.io/biocontainers/wget:1.20.1'

    publishDir "${params.db_dir}", mode: 'copy'

    output:
    path "dbcan/"

    script:
    def base = "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current"
    """
    mkdir -p dbcan
    cd dbcan

    wget_opts="-c --tries=0 --retry-connrefused --timeout=30 --read-timeout=60"

    wget \$wget_opts -O CAZy.dmnd              "${base}/CAZy.dmnd"              &
    wget \$wget_opts -O dbCAN.hmm              "${base}/dbCAN.hmm"              &
    wget \$wget_opts -O dbCAN-sub.hmm          "${base}/dbCAN_sub.hmm"          &
    wget \$wget_opts -O fam-substrate-mapping.tsv "${base}/fam-substrate-mapping.tsv" &
    wget \$wget_opts -O TCDB.dmnd              "${base}/TCDB.dmnd"              &
    wget \$wget_opts -O TF.hmm                 "${base}/TF.hmm"                 &
    wget \$wget_opts -O TF.dmnd                "${base}/TF.dmnd"                &
    wget \$wget_opts -O STP.hmm                "${base}/STP.hmm"                &
    wget \$wget_opts -O PUL.dmnd               "${base}/PUL.dmnd"               &
    wget \$wget_opts -O dbCAN-PUL.xlsx         "${base}/dbCAN-PUL.xlsx"         &
    wget \$wget_opts -O dbCAN-PUL.tar.gz       "${base}/dbCAN-PUL.tar.gz"       &
    wget \$wget_opts -O peptidase_db.dmnd      "${base}/peptidase_db.dmnd"      &
    wget \$wget_opts -O sulfatlas_db.dmnd      "${base}/sulfatlas_db.dmnd"      &
    wait

    gzip -dc dbCAN-PUL.tar.gz | tar -xf -
    """

    stub:
    """
    mkdir -p dbcan
    touch dbcan/CAZy.dmnd dbcan/dbCAN.hmm dbcan/dbCAN-sub.hmm
    touch dbcan/TCDB.dmnd dbcan/TF.hmm dbcan/TF.dmnd dbcan/STP.hmm dbcan/PUL.dmnd
    touch dbcan/dbCAN-PUL.xlsx dbcan/fam-substrate-mapping.tsv
    touch dbcan/peptidase_db.dmnd dbcan/sulfatlas_db.dmnd
    """
}

// ── 4. MEROPS pepunit database (peptidase annotation) ────────────────────────
// Downloads pepunit.lib from EBI MEROPS FTP and builds a DIAMOND index
process DOWNLOAD_MEROPS {
    tag "merops"
    label 'process_high'
    container 'quay.io/biocontainers/diamond:2.1.10--h43eeafb_0'

    publishDir "${params.db_dir}", mode: 'copy'

    output:
    path "merops/"

    script:
    """
    mkdir -p merops
    wget -q https://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib \\
        -O merops/merops_pepunit_raw.fasta
    # pepunit.lib has spaces in sequence lines — strip them before indexing
    sed '/^>/!s/ //g' merops/merops_pepunit_raw.fasta > merops/merops_pepunit.fasta
    rm merops/merops_pepunit_raw.fasta
    diamond makedb \\
        --in merops/merops_pepunit.fasta \\
        -d merops/merops_pepunit \\
        --threads ${task.cpus} \\
        --quiet
    """

    stub:
    """
    mkdir -p merops
    touch merops/merops_pepunit.fasta merops/merops_pepunit.dmnd
    """
}

// ── 5. PHI-base (pathogen-host interaction database) ─────────────────────────
// Downloads phi-base_current.fas from GitHub and builds a DIAMOND index
process DOWNLOAD_PHI_BASE {
    tag "phi_base"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.1.10--h43eeafb_0'

    publishDir "${params.db_dir}", mode: 'copy'

    output:
    path "phi_base/"

    script:
    """
    mkdir -p phi_base
    wget -q https://github.com/PHI-base/data/raw/master/releases/phi-base_current.fas \\
        -O phi_base/phi_base.fas
    diamond makedb \\
        --in phi_base/phi_base.fas \\
        -d phi_base/phi_base \\
        --threads ${task.cpus} \\
        --quiet
    """

    stub:
    """
    mkdir -p phi_base
    touch phi_base/phi_base.fas phi_base/phi_base.dmnd
    """
}

// ── 6. Gene Ontology OBO ──────────────────────────────────────────────────────
process DOWNLOAD_GO_OBO {
    tag "go_obo"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.11'

    publishDir "${params.db_dir}", mode: 'copy'

    output:
    path "go-basic.obo"

    script:
    """
    wget -q http://purl.obolibrary.org/obo/go/go-basic.obo -O go-basic.obo
    """

    stub:
    """
    touch go-basic.obo
    """
}

// ── 7. gene2product curated name database ─────────────────────────────────
// nextgenusfs/gene2product: 34k gene-name → NCBI-compliant product mappings.
// Used by ANNOTATE_FUNCTIONAL to normalise product names.
process DOWNLOAD_GENE2PRODUCT {
    tag "gene2product"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.11'

    publishDir "${params.db_dir}", mode: 'copy'

    output:
    path "ncbi_cleaned_gene_products.txt"

    script:
    """
    wget -q https://raw.githubusercontent.com/nextgenusfs/gene2product/master/ncbi_cleaned_gene_products.txt \\
        -O ncbi_cleaned_gene_products.txt
    """

    stub:
    """
    printf '#version 1.0\\n#Name\\tDescription\\nACT1\\tActin-1\\nTUB2\\tTubulin beta chain\\n' \\
        > ncbi_cleaned_gene_products.txt
    """
}

// ── 8. InterProScan 6 member databases ───────────────────────────────────────
// Downloads all member databases that require local data files (has_data = true
// in applications.config): AntiFam, CATH, CDD, HAMAP, NCBIFAM, PANTHER, Pfam,
// PIRSF, PIRSR, PRINTS, PROSITE, SFLD, SMART, SUPERFAMILY (~50 GB total).
// The interpro metadata archive is also downloaded (GO terms, pathways, entries).
// Extracted under interproscan/ to match the path IPS6 expects when
// --IPS6_databases_path <db_dir>/interproscan is given.
process DOWNLOAD_INTERPROSCAN {
    tag "interproscan"
    label 'process_long'
    container 'quay.io/biocontainers/python:3.11'

    publishDir "${params.db_dir}", mode: 'copy'

    output:
    path "interproscan/"

    script:
    """
    python3 download_interproscan.py
    """

    stub:
    """
    mkdir -p interproscan/interpro/109.0
    touch interproscan/interpro/109.0/databases.json
    for db in antifam/8.0 cath/4.3.0 cdd/3.21 hamap/2026_01 ncbifam/19.0 \\
              panther/19.0 pfam/38.2 pirsf/3.10 pirsr/2025_05 \\
              prints/42.0 prosite/2026_01 sfld/4 smart/9.0 superfamily/1.75; do
        mkdir -p interproscan/\${db}
        touch interproscan/\${db}/.done
    done
    """
}

// ── DOWNLOAD_DATABASES workflow ───────────────────────────────────────────────
workflow DOWNLOAD_DATABASES {

    if (!params.db_dir) {
        log.error "SETUP requires '--db_dir <path>' — specify where to save the databases."
        System.exit(1)
    }

    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║        bagRNA database setup                         ║
    ╚══════════════════════════════════════════════════════╝
    Download target   : ${params.db_dir}
    Skip eggNOG       : ${params.skip_eggnog_db}
    Skip Rfam         : ${params.skip_rfam_db}
    Skip dbCAN        : ${params.skip_dbcan_db}
    Skip MEROPS       : ${params.skip_merops_db}
    Skip PHI-base     : ${params.skip_phi_base_db}
    Skip GO OBO       : ${params.skip_go_obo_db}
    Skip gene2product : ${params.skip_gene2product_db}
    Skip InterProScan : ${params.skip_interproscan_db}
    """

    if (!params.skip_eggnog_db)       DOWNLOAD_EGGNOG()
    if (!params.skip_rfam_db)         DOWNLOAD_RFAM()
    if (!params.skip_dbcan_db)        DOWNLOAD_DBCAN()
    if (!params.skip_merops_db)       DOWNLOAD_MEROPS()
    if (!params.skip_phi_base_db)     DOWNLOAD_PHI_BASE()
    if (!params.skip_go_obo_db)        DOWNLOAD_GO_OBO()
    if (!params.skip_gene2product_db)  DOWNLOAD_GENE2PRODUCT()
    if (!params.skip_interproscan_db)  DOWNLOAD_INTERPROSCAN()
}
