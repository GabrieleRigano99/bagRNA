process MERGE_ADDITIONAL_MODELS {
    tag "merge_additional_models"
    label 'process_medium'
    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1'

    publishDir "${params.outdir}/structural_annotation", mode: 'copy'

    input:
    path annevo_gff
    path barrnap_gff
    path trnascan_gff
    path noncoding_gff
    path coding_mikado_gff
    path helixer_gff
    path liftover_gff      // pass empty file if absent

    output:
    path "evm_merged.gff", emit: merged_gff

    shell:
    """
    # Mikado-first: strip codon/UTR features from renamed Mikado output and use as base.
    # ANNEVO, Helixer, and other sources fill in loci not covered by Mikado.
    grep -v '^\s*#' !{coding_mikado_gff} \\
        | awk 'BEGIN{FS=OFS="\\t"} \$3 !~ "codon" && \$3 !~ "UTR"' \\
        > evm_merged.gff

    # Append non-overlapping features from a secondary GFF into evm_merged.gff
    append_non_overlapping() {
        local src="\$1"
        [ -s "\$src" ] || return 0
        local ids
        ids=\$(bedtools intersect \\
            -a "\$src" \\
            -b evm_merged.gff \\
            -wa -s -v \\
        | awk 'BEGIN{FS=OFS="\\t"} \$3=="mRNA" || \$3=="tRNA" || \$3=="rRNA" || \$3=="gene" {
            n=split(\$9, attrs, ";");
            for (i=1; i<=n; i++) { if (attrs[i] ~ /^ID=/) { print substr(attrs[i], 4); break } }
        }')
        [ -z "\$ids" ] && return 0
        echo "\$ids" \\
        | grep -F -w -f - "\$src" \\
        | sed 's/Alias.*//g; s/Name.*\\t//g' \\
        | awk '\$0 !~ "#" && \$3 !~ "intron" && \$3 !~ "codon" && \$3 !~ "UTR"' \\
        | sed 's/transcript/mRNA/g; s/other_pred1/Helixer/g; s/pasa/Mikado_loci/g; s/AGAT/Barrnap/g' \\
        >> evm_merged.gff
    }

    # rRNA/tRNA come from dedicated, independently-validated single-source
    # callers (Barrnap / tRNAscan-SE) rather than speculative gene-model
    # candidates. Unlike ANNEVO/Helixer/BUSCO-recovery loci, their overlapping
    # a protein-coding call is normal biology (tRNAs commonly sit inside
    # introns/UTRs) rather than redundant evidence, so they are appended
    # unconditionally instead of being dropped by the overlap filter.
    # tRNAscan-SE's own "pseudogene" calls are excluded (not real tRNA genes).
    append_all() {
        local src="\$1"
        [ -s "\$src" ] || return 0
        awk 'BEGIN{FS=OFS="\\t"} \$3=="pseudogene"{
            n=split(\$9, attrs, ";");
            for (i=1; i<=n; i++) { if (attrs[i] ~ /^ID=/) { print substr(attrs[i], 4); break } }
        }' "\$src" > .pseudo_ids.txt
        grep -v '^\s*#' "\$src" \\
        | if [ -s .pseudo_ids.txt ]; then grep -v -w -F -f .pseudo_ids.txt; else cat; fi \\
        | sed 's/Alias.*//g; s/Name.*\\t//g' \\
        | awk '\$3 !~ "intron" && \$3 !~ "codon" && \$3 !~ "UTR"' \\
        | sed 's/transcript/mRNA/g; s/other_pred1/Helixer/g; s/pasa/Mikado_loci/g; s/AGAT/Barrnap/g' \\
        >> evm_merged.gff
    }

    append_non_overlapping !{annevo_gff}
    append_non_overlapping !{noncoding_gff}
    append_non_overlapping !{helixer_gff}
    append_all !{barrnap_gff}
    append_all !{trnascan_gff}

    append_non_overlapping !{liftover_gff}
    """

    stub:
    """
    touch evm_merged.gff
    """
}
