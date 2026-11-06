#!/bin/sh
# NCBI rejects EC_number on a "hypothetical protein" (SEQ_FEAT.BadProteinName). table2asn -c f
# renames some products (e.g. "protein") to "hypothetical protein" AFTER the .tbl is written, so
# the EC guard in gff3_to_tbl.py cannot see it. Usage: strip_ec_from_renamed.sh in.gbf in.tbl out.tbl
# Prints the number of CDS that lost their EC_number(s).
gbf=$1; tbl=$2; out=$3
awk 'function flush(){ if (blk && ec && hyp && lt != "") print lt }
     /^     [A-Za-z_]/ { flush(); blk = ($1 == "CDS"); ec = 0; hyp = 0; lt = "" }
     blk { if ($0 ~ /\/EC_number=/) ec = 1
           if ($0 ~ /\/product="hypothetical protein"/) hyp = 1
           if (match($0, /\/locus_tag="[^"]+"/)) lt = substr($0, RSTART + 12, RLENGTH - 13) }
     END { flush() }' "$gbf" | sort -u > bad_loci.txt
awk 'NR == FNR { bad[$1] = 1; next }
     /^[<>0-9]/ { drop = 0 }
     $1 == "locus_tag" && ($2 in bad) { drop = 1 }
     $1 == "EC_number" && drop { next }
     { print }' bad_loci.txt "$tbl" > "$out"
wc -l < bad_loci.txt
