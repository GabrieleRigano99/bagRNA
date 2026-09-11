process FILTER_ISOFORMS {
    tag "filter_isoforms"
    label 'process_low'
    container 'ubuntu:22.04'

    publishDir "${params.outdir}/mikado", mode: 'copy'

    input:
    path pick_gff

    output:
    path "isoform_filtered.gff3", emit: filtered_gff

    script:
    // Keep primary=True transcripts always.
    // Keep primary=False only when ccode=j (genuine AS isoform with different junction).
    // Drop all other primary=False (redundant coverage, retained-intron artefacts).
    """
    awk '
    BEGIN { skip_id = "" }

    /^#/ { print; next }

    \$3 == "mRNA" {
        if (index(\$9, "primary=False") > 0 && index(\$9, "ccode=j") == 0) {
            n = split(\$9, a, ";")
            for (i = 1; i <= n; i++) {
                if (index(a[i], "ID=") == 1) { skip_id = substr(a[i], 4); break }
            }
        } else {
            skip_id = ""
            print
        }
        next
    }

    \$3 == "gene" || \$3 == "superlocus" || \$3 == "ncRNA_gene" {
        skip_id = ""; print; next
    }

    skip_id != "" {
        n = split(\$9, a, ";")
        for (i = 1; i <= n; i++) {
            if (index(a[i], "Parent=") == 1 && substr(a[i], 8) == skip_id) next
        }
        print; next
    }

    { print }
    ' ${pick_gff} > isoform_filtered.gff3
    """

    stub:
    """
    touch isoform_filtered.gff3
    """
}
