#!/bin/bash

# Default settings
stringtie_gtf=""
trinity_transcripts=""
gff=""
fasta=""
locus_tag=""
interproscan_sh_path=""
signalp_models_path=""
EggNogdb_path=""
C1_R1=""
C1_R2=""
C2_R1=""
C2_R2=""
C3_R1=""
C3_R2=""
T1_R1=""
T1_R2=""
T2_R1=""
T2_R2=""
T3_R1=""
T3_R2=""
intermediate_files_dir="intermediate_files"
annotation_files_dir="annotation_files"

# Function to print the correct usage of the script
usage() {
    echo "Usage: $0 --stringtie_gtf <file> --trinity_transcripts <file> --gff <file> --fasta <file> --locus_tag <tag> --interproscan_sh_path <path> --signalp_models_path <path> --EggNogdb_path <path> --C1_R1 <file> --C1_R2 <file> --C2_R1 <file> --C2_R2 <file> --C3_R1 <file> --C3_R2 <file> --T1_R1 <file> --T1_R2 <file> --T2_R1 <file> --T2_R2 <file> --T3_R1 <file> --T3_R2 <file>"
    echo "    --stringtie_gtf <file>: Specify the input GTF file from Stringtie"
    echo "    --trinity_transcripts <file>: Specify the input FASTA file from Trinity"
    echo "    --gff <file>: Specify the input GFF file from prediction"
    echo "    --fasta <file>: Specify the input FASTA file for the genome"
    echo "    --locus_tag <tag>: Specify the locus tag"
    echo "    --interproscan_sh_path <path>: Specify the path to interproscan.sh"
    echo "    --signalp_models_path <path>: Specify the path to the SignalP models directory"
    echo "    --EggNogdb_path <path>: Specify the path to the EggNOG database"
    echo "    --C1_R1 <file>: Path to R1 file for sample C1"
    echo "    --C1_R2 <file>: Path to R2 file for sample C1"
    echo "    --C2_R1 <file>: Path to R1 file for sample C2"
    echo "    --C2_R2 <file>: Path to R2 file for sample C2"
    echo "    --C3_R1 <file>: Path to R1 file for sample C3"
    echo "    --C3_R2 <file>: Path to R2 file for sample C3"
    echo "    --T1_R1 <file>: Path to R1 file for sample T1"
    echo "    --T1_R2 <file>: Path to R2 file for sample T1"
    echo "    --T2_R1 <file>: Path to R1 file for sample T2"
    echo "    --T2_R2 <file>: Path to R2 file for sample T2"
    echo "    --T3_R1 <file>: Path to R1 file for sample T3"
    echo "    --T3_R2 <file>: Path to R2 file for sample T3"
    echo "    --intermediate_files_dir <directory>: Specify the name of the output folder for intermediate files (default: intermediate_files)"
    echo "    --annotation_files_dir <directory>: Specify the name of the output folder for annotation files (default: annotation_files)"
    exit 1
}

# Parsing command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --stringtie_gtf)
            stringtie_gtf="$2"
            shift
            shift
            ;;
        --trinity_transcripts)
            trinity_transcripts="$2"
            shift
            shift
            ;;
        --gff)
            gff="$2"
            shift
            shift
            ;;
        --fasta)
            fasta="$2"
            shift
            shift
            ;;
        --locus_tag)
            locus_tag="$2"
            shift
            shift
            ;;
        --interproscan_sh_path)
            interproscan_sh_path="$2"
            shift
            shift
            ;;
        --signalp_models_path)
            signalp_models_path="$2"
            shift
            shift
            ;;
        --EggNogdb_path)
            EggNogdb_path="$2"
            shift
            shift
            ;;
        --intermediate_files_dir)
            intermediate_files_dir="$2"
            shift
            shift
            ;;
        --annotation_files_dir)
            annotation_files_dir="$2"
            shift
            shift
            ;;
           --C1_R1)
            C1_R1="$2"
            shift
            shift
            ;;
        --C1_R2)
            C1_R2="$2"
            shift
            shift
            ;;
        --C2_R1)
            C2_R1="$2"
            shift
            shift
            ;;
        --C2_R2)
            C2_R2="$2"
            shift
            shift
            ;;
        --C3_R1)
            C3_R1="$2"
            shift
            shift
            ;;
        --C3_R2)
            C3_R2="$2"
            shift
            shift
            ;;
        --T1_R1)
            T1_R1="$2"
            shift
            shift
            ;;
        --T1_R2)
            T1_R2="$2"
            shift
            shift
            ;;
        --T2_R1)
            T2_R1="$2"
            shift
            shift
            ;;
        --T2_R2)
            T2_R2="$2"
            shift
            shift
            ;;
        --T3_R1)
            T3_R1="$2"
            shift
            shift
            ;;
        --T3_R2)
            T3_R2="$2"
            shift
            shift
            ;;
        *)
            usage
            ;;
    esac
done

# Check if all necessary inputs are provided
if [ -z "$stringtie_gtf" ] || [ -z "$trinity_transcripts" ] || [ -z "$gff" ] || [ -z "$fasta" ] || [ -z "$locus_tag" ] || [ -z "$C1_R1" ] || [ -z "$C1_R2" ] || [ -z "$C2_R1" ] || [ -z "$C2_R2" ] || [ -z "$C3_R1" ] || [ -z "$C3_R2" ] || [ -z "$T1_R1" ] || [ -z "$T1_R2" ] || [ -z "$T2_R1" ] || [ -z "$T2_R2" ] || [ -z "$T3_R1" ] || [ -z "$T3_R2" ] || [ -z "$interproscan_sh_path" ] || [ -z "$signalp_models_path" ] || [ -z "$EggNogdb_path" ]; then    
    usage
fi

# Create folder for intermediate files and annotation files  
mkdir -p "$intermediate_files_dir"
mkdir -p "$annotation_files_dir"

message="Ma cumminamu un maciellu?"
sleep 5 && for ((i=0; i<${#message}; i++)); do echo -n -e "\e[1;31m${message:$i:1}\e[0m"; sleep 0.1; done && sleep 3

# rRNA and tRNA annotation

barrnap "$fasta" --kingdom euk --threads 12 > barrnap.gff
tRNAscan-SE "$fasta" --thread 12 -j tRNA.gff -Q

# Generate, convert, complement and clean annotations
funannotate stringtie2gff3 -i "$stringtie_gtf" > Stringtie.gff
gmap_build -d Trinity -D . -k 13 "$fasta"
gmap -t 12 -n 0 -D ./ -d Trinity -f gff3_gene "$trinity_transcripts" > Trinity.gff
agat_sp_complement_annotations.pl --ref Stringtie.gff --add Trinity.gff -o complement_stringtie_trinity.gff
bedtools intersect -a "$gff" -b complement_stringtie_trinity.gff -v -s > tmp1; cat tmp1 complement_stringtie_trinity.gff > intersect_all.gff; rm tmp1
agat_sp_fix_overlaping_genes.pl -f intersect_all.gff -o olp_intersect.gff
agat_sp_fix_features_locations_duplicated.pl -f olp_intersect.gff -o locdup_olp_intersect.gff
funannotate gff-rename -g locdup_olp_intersect.gff -f "$fasta" -o rename1.gff -l "$locus_tag"

# Transuite evaluation of gene models 

agat_convert_sp_gff2gtf.pl --gff structuralv1.gff -o structuralv1.gtf
agat_sp_extract_sequences.pl --mrna -g structuralv1.gtf -f "$fasta" -o transuite_mrna.fa
python ~/tools/TranSuite/transuite.py Auto --gtf structuralv1.gtf --fasta transuite_mrna.fa --outpath ./ --outname transuite --cds 30 --iter 5 --pep 100 --ptc 70
agat_convert_sp_gxf2gxf.pl -g transuite_TranSuite_output/transuite_transfix/transuite_transfix.gtf -o transuite_transfix.gff
awk 'BEGIN{FS=OFS="\t"}{if ($3=="RNA") gsub (/RNA/, "ncRNA", $0); print }' transuite_transfix.gff > transuite.gff
agat_sp_fix_cds_phases.pl --gff transuite.gff -f "$fasta" -o cds_fixed.gff
funannotate gff-rename -g cds_fixed.gff -f "$fasta" -o rename2.gff -l "$locus_tag"

# remove ncRNAs below 0.1 TPM
sed 's/Alias=/\t/g' rename2.gff | cut -f1-9 > structuralv2.gff
agat_convert_sp_gff2gtf.pl --gff structuralv2.gff -o structuralv2.gtf
agat_sp_extract_sequences.pl --mrna --gff structuralv2.gff --fasta "$fasta" -o transcripts1.fa
~/tools/CPC2_standalone-1.0.1/bin/CPC2.py -i transcripts1.fa -o CPC2_result
awk 'BEGIN{FS=OFS="\t"}{print $1,$8}' CPC2_result.txt | sed 's/-T[0-9]*//g' | bedtools groupby -g 1 -c 2 -o collapse > genes_list_CPC2.tsv
grep -w "coding" genes_list_CPC2.tsv | cut -f1 > geni_codingv1.tsv
grep -v -f geni_codingv1.tsv genes_list_CPC2.tsv | cut -f1 > geni_noncodingv1.tsv
salmon index --transcripts transcripts1.fa -p 12 --index transcripts_index
salmon quant --seqBias --posBias --gcBias -g structuralv2.gtf -p 12 -l A -i transcripts_index -1 "$C1_R1" -2 "$C1_R2" --output C1
salmon quant --seqBias --posBias --gcBias -g structuralv2.gtf -p 12 -l A -i transcripts_index -1 "$C2_R1" -2 "$C2_R2" --output C2
salmon quant --seqBias --posBias --gcBias -g structuralv2.gtf -p 12 -l A -i transcripts_index -1 "$C3_R1" -2 "$C3_R2" --output C3
salmon quant --seqBias --posBias --gcBias -g structuralv2.gtf -p 12 -l A -i transcripts_index -1 "$T1_R1" -2 "$T1_R2" --output T1
salmon quant --seqBias --posBias --gcBias -g structuralv2.gtf -p 12 -l A -i transcripts_index -1 "$T2_R1" -2 "$T2_R2" --output T2
salmon quant --seqBias --posBias --gcBias -g structuralv2.gtf -p 12 -l A -i transcripts_index -1 "$T3_R1" -2 "$T3_R2" --output T3
salmon quantmerge --genes --quants C1 C2 C3 T1 T2 T3 --names C1 C2 C3 T1 T2 T3 -o quant.tsv
awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\n",$1,$2,$3,$4,$5,$6,$7,($2+$3+$4)/3,($5+$6+$7)/3}' quant.tsv | awk 'BEGIN{FS=OFS="\t"}{if($8>=0.1 || $9>=0.1) print}' > filtered_quant_genes_0.1TPM.tsv
grep -f geni_noncodingv1.tsv filtered_quant_genes_0.1TPM.tsv | cut -f1 > noncoding_over_cutoff.tsv
grep -v -f noncoding_over_cutoff.tsv geni_noncodingv1.tsv > ncs_to-remove.tsv
grep -v -f ncs_to-remove.tsv structuralv2.gff > structuralv3.gff

# complement tRNAs and rRNAs to structural and rename final annotation

agat_sp_merge_annotations.pl -f barrnap.gff -f tRNA.gff -f structuralv3.gff -o structuralv4.gff
funannotate gff-rename -g structuralv4.gff -f "$fasta" -o rename3.gff -l "$locus_tag"
sed 's/;Name=SS[0-9]*_[0-9]*.tRNA[0-9]*/;Name=tRNA/g' rename3.gff | sed 's/;product=None//g' | sed 's/Alias/\t/g' | cut -f1-9 > final_structural.gff3

# update coding and noncoding files

agat_sp_extract_sequences.pl --mrna --gff final_structural.gff3 --fasta "$fasta" -o transcripts_final_structural.fa
~/tools/CPC2_standalone-1.0.1/bin/CPC2.py -i transcripts_final_structural.fa -o CPC2_result2
awk 'BEGIN{FS=OFS="\t"}{print $1,$8}' CPC2_result2.txt | sed 's/-T[0-9]*//g' | bedtools groupby -g 1 -c 2 -o collapse > genes_list_CPC2_2.tsv
grep -w "coding" genes_list_CPC2_2.tsv | cut -f1 > geni_coding_final_structural.txt
grep -v -f geni_coding_final_structural.txt genes_list_CPC2_2.tsv | cut -f1 > geni_noncodingv2.tsv
awk -F "\t" '{if($3=="tRNA" || $3=="rRNA") print $9}' final_structural.gff3 | sed 's/ID=//g' | sed 's/-/\t/g' | cut -f1 > tRNA_rRNA_list.tsv
grep -v -f tRNA_rRNA_list.tsv geni_noncodingv2.tsv > geni_noncoding_final_structural.txt

mv *.gff *CPC2* *.tsv *.fa *.gtf C1 C2 C3 T1 T2 T3 Trinity transcripts_index transuite_TranSuite_output *report* *.log "$intermediate_files_dir"
# Functional annotation with IPRscan, EggNog, SignalP6.0, phobius

funannotate gff2prot -g final_structural.gff3 -f "$fasta" --no_stop > proteins.fasta
fasta2tab.awk proteins.fasta | awk 'BEGIN{FS="\t"; OFS="\n"}{if($2!~"*") print ">"$1,$2}' > proteins.no_inner_stop.fasta
funannotate iprscan -i proteins.no_inner_stop.fasta -m local -o IPRscan.XML -c 10 --iprscan_path "$interproscan_sh_path"
signalp6 --fastafile proteins.no_inner_stop.fasta --organism eukarya --output_dir signalp_op --format txt --mode slow-sequential -md "$signalp_models_path" -tt 10 -wp 10 
~/miniconda3/envs/funannotate/bin/emapper.py -i proteins.fasta --output eggnog --data_dir "$EggNogdb_path" --cpu 12 --override
~/tools/phobius101_linux/phobius/phobius.pl proteins.no_inner_stop.fasta -short > phobius.txt 
echo lanciare BlastKoala su https://www.kegg.jp/blastkoala/ e poi lanciare il successivo script per finire la funzionale

mv *.fasta *eggnog* *.XML *signalp* *.gff3 *.txt "$annotation_files_dir"

message="Ma cumminamu un maciellu?"
sleep 5 && for ((i=0; i<${#message}; i++)); do echo -n -e "\e[1;31m${message:$i:1}\e[0m"; sleep 0.1; done && sleep 3