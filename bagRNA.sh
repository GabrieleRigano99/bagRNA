#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

# Save original stdout/stderr before redirecting
exec 3>&1
exec 4>&2

# Trap for general errors
trap 'echo -e "\e[31m[ERROR]\e[0m An unexpected error occurred. Exiting." >&2' ERR

# Create or clear log files at start
: > pipeline.log
: > error.log

# Log stdout to pipeline.log
exec > >(ts '[%Y-%m-%d %H:%M:%S]' | tee -a pipeline.log)

# Log stderr to both terminal and error.log
exec > >(ts '[%Y-%m-%d %H:%M:%S]' | tee -a error.log >&2)

# Default settings
threads="12"
busco_lineage=""
input_fasta=""
locus_tag="bagRNA"
Helixer_gff=""
RAM_limit_Trinity="45G"
limitBAMsortRAM="2162824741"
orientation=""
strandedness=""
max_intron_length=""
prot_evidence=""
codon_table="1"
species=""
scoring=""
mikado_config=""
lifted_annotation=""
STAR_manifest=()
databases=""
strain="strain"
submission_template=""
no_functional_anno=false
jaccard_clip=true
no_antismash=false
no_effectorp3=false
no_tmbed=false
GeneMark_PATH=""
phobius_PATH=""
Helixer_lineage=""
max_gene_length=""

# === Fixed usage function ===
usage() {
    # Restore terminal output
    exec 1>&3 >&4
 
    cat >&3 << 'EOF'

        # =============================================== #
        #   _                                             # 
        #  | |                  ____    _   _      _      #
        #  | |_    ____  __ _  |  _ \  | \ | |    /_\     #
        #  |  _ \ / _` |/ _` | | |_) | |  \| |   //_\\    #
        #  | |_) | (_| | (_| | |  _ <  | |\  |  / ___ \   #
        #  |____/ \__,_|\__, | |_| \_\ |_| \_| /_/   \_\  #
        #               |___/                             #
        #                                                 #
        # =============================================== #

EOF

message="
"

sleep 1 && for ((i=0; i<${#message}; i++)); do echo -n -e "\e[1;31m${message:$i:1}\e[0m"; sleep 0.1; done && sleep 1

    echo "your end-to-end noncoding-aware eukaryotic genome annotation pipeline " >&3
    echo "" >&3
    echo "bagRNA version 1.0.0">&3
    echo "" >&3
    echo "Usage:         bagRNA <arguments>" >&3
    echo "" >&3
    echo "Description:   The bagRNA annotation pipeline takes as input a genome, proteins and RNA-seq short reads (from the STAR tsv manifest file)
               to generate a full annotation of protein-coding and non-coding genes together with their functions.
               The user is required to provide configuration files specified in the mikado documentation and a submission template file 
               (examples are present in the bagRNA GitHub page https://github.com/BCGL-Unime/bagRNA)
               Written and developed by Gabriele Rigano - Bioinformatics and Computational Genomics LAB - University of Messina - gabrielerigano99@gmail.com ">&3
    echo "" >&3
    echo "  -h, --help             Display this help message and exit" >&3
    echo "" >&3
    echo "Mandatory inputs:" >&3
    echo "" >&3
    echo "  --input_fasta <file>         Genome FASTA file (better if softmasked)" >&3
    echo "  --prot_evidence <file>       Protein evidence in FASTA format" >&3
    echo "  --busco_lineage <string>     BUSCO lineage (e.g. sordariomycetes)" >&3
    echo "  --STAR_manifest <tsv file/s> STAR manifest TSV files for STAR_manifest, an example is provided in the bagRNA github repository (https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf)  " >&3
    echo "  --mikado_config <file>       Specify the mikado configuration file, an example is provided in the bagRNA github repository (https://mikado.readthedocs.io/en/stable/Tutorial/)" >&3
    echo "  --scoring <string>           Scoring file for mikado (e.g. scerevisiae.yaml) (https://mikado.readthedocs.io/en/stable/Tutorial/Scoring_tutorial/)" >&3
    echo "  --species <string>           Species name, must be written in quotes with underscores instead of spaces (e.g. \"Arabidopsis_thaliana\")" >&3
    echo "  --submission_template <file> Sbt submission template .sbt, an example is provided in the bagRNA github repository" >&3
    echo "" >&3
    echo "Choose one flag among:" >&3
    echo ""
    echo "  --Helixer_lineage <string>   Lineage for helixer prediction (requires GPU and NVIDIA Container Toolkit configuration) (choose among: fungi, land_plant, vertebrate, invertebrate)" >&3
    echo "  --Helixer_gff <file>         Input GFF file from precomputed Helixer run (recommended) (https://www.plabipd.de/helixer_main.html)" >&3
    echo ""
    echo "Additional inputs:" >&3
    echo "" >&3
    echo "  --lifted_annotation <file>   Input GFF lifted annotation from liftoff/lifton" >&3
    echo "  --GeneMark_PATH <path>       Path to GeneMarK-ET/ETP (where the gmes_petap.pl executable is)" >&3
    echo "  --phobius_PATH <path>        Path to phobius executable (where the phobius.pl executable is)" >&3
    echo "  --databases <path>           Path to databases for functional annotation (databases can be downloaded with download_databases_bagRNA.sh script present in the bagRNA github repository)" >&3
    echo "" >&3
    echo "Performance and Miscellaneous options:" >&3
    echo "" >&3
    echo "  --threads, -t <int>          Number of threads (default: 20)" >&3
    echo "  --jaccard_clip               Use this flag if you are expecting high gene density with UTR overlap (recommended for fungi) (default: uses --jaccard_clip)" >&3
    echo "  --max_gene_length <int>      Maximum length a gene model can be (default: 30000 bp)" >&3
    echo "  --RAM_limit_Trinity <int>    RAM limit for Trinity (default: 45G)" >&3
    echo "  --limitBAMsortRAM <int>      STAR BAM sort RAM limit" >&3
    echo "  --orientation <string>       Reads orientation (e.g. FR, RF)" >&3
    echo "  --strandedness <string>      Strandedness type (e.g. secondstrand) (https://chipster.csc.fi/manual/library-type-summary.html)" >&3
    echo "  --max_intron_length <int>    Max intron length (default: 3000)" >&3
    echo "  --codon_table <int>          Codon table (default: 1)" >&3
    echo "  --strain <string>            Strain/Isolate name (default: strain)" >&3
    echo "  --locus_tag <tag>            Locus tag (default: bagRNA)" >&3
    echo "  --no_functional_anno         Use this flag if you do NOT wish to run functional annotation (default: runs functional annotation)" >&3
    echo "  --no_antismash               Use this flag if you do NOT wish to run AntiSMASH (default: runs AntiSMASH)" >&3
    echo "  --no_effectorp3              Use this flag if you do NOT wish to run EffectorP-3.0 (default: runs EffectorP-3.0)" >&3    
    echo "  --no_tmbed                   Use this flag if you do NOT wish to run Tmbed (default: runs Tmbed)" >&3
    echo "" >&3
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -h|--help)
            usage
            ;;
        --busco_lineage)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --busco_lineage" && usage
            busco_lineage="$2"
            shift; shift
            ;;
        --Helixer_gff)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --Helixer_gff" && usage
            Helixer_gff="$2"
            shift; shift
            ;;
        --prot_evidence)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --prot_evidence" && usage
            prot_evidence="$2"
            shift; shift
            ;;
        --input_fasta)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --input_fasta" && usage
            input_fasta="$2"
            shift; shift
            ;;
        --RAM_limit_Trinity)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --RAM_limit_Trinity" && usage
            RAM_limit_Trinity="$2"
            shift; shift
            ;;
        --limitBAMsortRAM)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --limitBAMsortRAM" && usage
            limitBAMsortRAM="$2"
            shift; shift
            ;;
        --locus_tag)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --locus_tag" && usage
            locus_tag="$2"
            shift; shift
            ;;
        --orientation)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --orientation" && usage
            orientation="$2"
            shift; shift
            ;;
        --strandedness)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --strandedness" && usage
            strandedness="$2"
            shift; shift
            ;;
        --max_intron_length)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --max_intron_length" && usage
            max_intron_length="$2"
            shift; shift
            ;;
        --STAR_manifest)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                STAR_manifest+=("$1")
                shift
            done
            ;;
        -t|--threads)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --threads" && usage
            threads="$2"
            shift; shift
            ;;
        --species)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --species" && usage
            species="$2"
            shift; shift
            ;;
        --mikado_config)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --mikado_config" && usage
            mikado_config="$2"
            shift; shift
            ;;
        --codon_table)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --codon_table" && usage
            codon_table="$2"
            shift; shift
            ;;
        --scoring)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --scoring" && usage
            scoring="$2"
            shift; shift
            ;;
        --databases)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --databases" && usage
            databases="$2"
            shift; shift
            ;;
        --strain)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --strain" && usage
            strain="$2"
            shift; shift
            ;;
        --submission_template)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --submission_template" && usage
            submission_template="$2"
            shift; shift
            ;;
        --Helixer_lineage)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --Helixer_lineage" && usage
            Helixer_lineage="$2"
            shift; shift
            ;;
        --max_gene_length)
            [[ -z "${2:-}" || "$2" == --* ]] && echo "Missing argument for --max_gene_length" && usage
            max_gene_length="$2"
            shift; shift
            ;;    
        --GeneMark_PATH)
            [[ -z "${2:-}" || "$2" == --* ]] && usage
            echo "GeneMark executable provided, bagRNA will use it in the gene prediction step"
            GeneMark_PATH="$2"
            shift; shift
            ;;
        --phobius_PATH)
            [[ -z "${2:-}" || "$2" == --* ]] && usage
            echo "phobius executable provided, bagRNA will use it in the functional annotation step"
            phobius_PATH="$2"
            shift; shift
            ;;            
        --jaccard_clip)
            echo "--jaccard_clip flag used, expecting high gene density with UTR overlap - recommended for fungi"
            jaccard_clip=true
            shift
            ;;
        --lifted_annotation)
            [[ -z "${2:-}" || "$2" == --* ]] && usage
            echo "lifted_annotation provided, bagRNA will use it in the gene prediction step"
            lifted_annotation="$2"
            shift; shift
            ;;
        --no_functional_anno)
            echo "--no_functional_anno flag used, bagRNA will skip functional annotation"
            no_functional_anno=true
            shift
            ;;
        --no_effectorp3)
            echo "--no_effectorp3 flag used, bagRNA will skip EffectorP-3.0 annotation"
            no_effectorp3=true
            shift
            ;;            
        --no_antismash)
            echo "--no_antismash flag used, bagRNA will skip AntiSMASH annotation"
            no_antismash=true
            shift
            ;;
        --no_tmbed)
            echo "--no_tmbed flag used, bagRNA will skip Tmbed annotation"
            no_tmbed=true
            shift
            ;;        
        *)
            echo "Error: Unrecognized argument '$1'"
            usage
            ;;
    esac
done

# === Validate required inputs ===
missing_inputs=()
required_vars=("busco_lineage" "input_fasta" "scoring" "mikado_config" "species" "prot_evidence" "submission_template")

for var in "${required_vars[@]}"; do
    # This checks if variable is either unset or empty
    if [ -z "${!var+x}" ] || [ -z "${!var}" ]; then
        missing_inputs+=("$var")
    fi
done

# STAR_manifest as array
if [ -z "${STAR_manifest+x}" ] || [ ${#STAR_manifest[@]} -eq 0 ]; then
    missing_inputs+=("STAR_manifest")
fi

if [ ${#missing_inputs[@]} -ne 0 ]; then
    echo -e "\n\e[31m[ERROR]\e[0m The following required inputs are missing:"
    for input in "${missing_inputs[@]}"; do
        echo "  - $input"
    done
    echo ""
    usage
fi


# Utility functions for colored output
print_red() {
    echo -e "\e[31m$1\e[0m"
}

print_green() {
    echo -e "\e[32m$1\e[0m"
}

# === Function to check if folders exist ===
check_folders_exist() {
    if [ -d tmp_dir ]  && [ -d working_dir ]; then
        print_green "Folders tmp_dir, working_dir and agat_config found. Skipping creation."
        return 0
    else
        print_red "Folders tmp_dir, working_dir and agat_config not found. Creating..."
        return 1
    fi
}

# === Function to check if BUSCO results exist ===
check_busco_results() {
    if [ -f "busco_genome/short_summary.specific."$busco_lineage"_odb12.busco_genome.txt" ] || [ -f "working_dir/busco_genome/short_summary.specific."$busco_lineage"_odb12.busco_genome.txt" ] ; then
        print_green "BUSCO results found. Skipping BUSCO run."
        return 0
    else
        print_red "BUSCO results not found. Running BUSCO..."
        return 1
    fi
}

# === Function to check if reformatted miniprot GFF exists ===
check_miniprot_results() {
    if [ -f "tmp_dir/busco_genome.gff" ]; then
        print_green "reformatted_miniprot.gff found. Skipping reformatting."
        return 0
    else
        print_red "reformatted_miniprot.gff not found. Running reformatting..."
        return 1
    fi
}

# === Function to check if fixed miniprot GFF exists ===
check_miniprot_fix() {
    if [ -f "working_dir/busco.gff" ]; then
        print_green "Fixed miniprot GFF found. Skipping agat_fix_cds."
        return 0
    else
        print_red "Fixed miniprot GFF not found. Running agat_fix_cds..."
        return 1
    fi
}

# === Function to check if Helixer GTF exists ===
check_helixer_gtf() {
    if [ -f "working_dir/helixer.gtf" ]; then
        print_green "Helixer GTF found. Skipping GFF to GTF conversion."
        return 0
    else
        print_red "Helixer GTF not found. Running GFF to GTF conversion..."
        return 1
    fi
}

# === Function to check if STAR genome index exists ===
check_star_genome_index() {
    if [ -d "tmp_dir/STAR_index_files" ]; then
        print_green "STAR genome index found. Skipping genomeGenerate."
        return 0
    else
        print_red "STAR genome index not found. Running genomeGenerate..."
        return 1
    fi
}

# === Function to check if all STAR BAMs exist ===
check_star_mapping() {
    local missing=0
    for i in "${STAR_manifest[@]}"; do
        name=$(basename "$i")
        if [ ! -f "${name}_Aligned.sortedByCoord.out.bam" ] && [ ! -f "tmp_dir/${name}_Aligned.sortedByCoord.out.bam" ]; then
            print_red "STAR mapping for condition $i not found. Running STAR mapping..."
            missing=1
        fi
    done
    if [ "$missing" -eq 0 ]; then
        print_green "All STAR mapping BAM files found. Skipping STAR mapping."
        return 0
    else
        return 1
    fi
}

# === Function to check if Portcullis filtered BAMs exist ===
check_portcullis_filtering() {
    local missing=0
    for i in "${STAR_manifest[@]}"; do
        name=$(basename "$i")
        if [ ! -f "Portcullis_"${name}"/portcullis.filtered.bam" ] && [ ! -f "working_dir/Portcullis_"${name}"/portcullis.filtered.bam" ]; then
            print_red "Portcullis filtered BAM for $name not found. Running Portcullis filtering..."
            missing=1
        fi
    done
    if [ "$missing" -eq 0 ]; then
        print_green "All Portcullis filtered BAMs found. Skipping Portcullis filtering."
        return 0
    else
        return 1
    fi
}

# === Function to check if merged Portcullis junctions BED exists ===
check_portcullis_merged() {
    if [ -f "working_dir/all_splice_junctions.bed" ]; then
        print_green "Merged Portcullis junctions found. Skipping merge."
        return 0
    else
        print_red "Merged Portcullis junctions not found. Running merge..."
        return 1
    fi
}

# === Function to check if merged BAM exists ===
check_bam_merge() {
    if [ -f "working_dir/merged_mappings.bam" ]; then
        print_green "Merged BAM file found. Skipping BAM merge."
        return 0
    else
        print_red "Merged BAM file not found. Running BAM merge..."
        return 1
    fi
}

# === Function to check if Aletsch assembly GTF exists ===
check_aletsch_assembly() {
    if [ -f "working_dir/aletsch.gtf" ]; then
        print_green "Aletsch transcript assembly found. Skipping Aletsch."
        return 0
    else
        print_red "Aletsch transcript assembly not found. Running Aletsch..."
        return 1
    fi
}

# === Function to check if StringTie assembly GTF exists ===
check_stringtie_assembly() {
    if [ -f "working_dir/stringtie.gtf" ]; then
        print_green "StringTie assembly found. Skipping StringTie assembly."
        return 0
    else
        print_red "StringTie assembly not found. Running StringTie assembly..."
        return 1
    fi
}

# === Function to check if reformatted Aletsch GTF exists ===
check_aletsch_merged() {
    if [ -f "working_dir/aletsch.gtf" ]; then
        print_green "Reformatted Aletsch GTF found. Skipping reformatting."
        return 0
    else
        print_red "Reformatted Aletsch GTF not found. Running reformatting..."
        return 1
    fi
}

# === Function to check if Trinity genome-guided assembly exists ===
check_trinity_assembly() {
    if [ -f "tmp_dir/op_trinity.Trinity-GG.fasta" ]; then
        print_green "Trinity assembly found. Skipping Trinity run."
        return 0
    else
        print_red "Trinity assembly not found. Running Trinity assembly..."
        return 1
    fi
}

# === Function to check if seqcleaned Trinity transcripts exist ===
check_seqcleaned_transcripts() {
    if [[ -f "op_trinity.Trinity-GG.fasta.clean" || -f "tmp_dir/op_trinity.Trinity-GG.fasta.clean" ]]; then
        print_green "Seqcleaned Trinity transcripts found. Skipping seqclean."
        return 0
    else
        print_red "Seqcleaned Trinity transcripts not found. Running seqclean..."
        return 1
    fi
}

# === Function to check if GMAP genome index exists ===
check_gmap_index() {
    if [ -d "tmp_dir/gmap_genome_index" ] || [ -d "gmap_genome_index" ]; then
        print_green "GMAP genome index found. Skipping gmap_index."
        return 0
    else
        print_red "GMAP genome index not found. Running gmap_index..."
        return 1
    fi
}

# === Function to check if GMAP GFF alignment exists ===
check_gmap_alignment() {
    if [ -f "working_dir/trinity.gff" ]; then
        print_green "Trinity alignment GFF found. Skipping gmap alignment."
        return 0
    else
        print_red "Trinity alignment GFF not found. Running gmap alignment..."
        return 1
    fi
}

# === Function to check if lifted annotation has been formatted ===
check_reformat_lifted_anno() {
    if [ -f "working_dir/liftover.gff" ]; then
        print_green "Reformatted lifted annotation GFF found. Skipping generation."
        return 0
    else
        print_red "Reformatted lifted annotation GFF not found. Reformatting..."
        return 1
    fi    
}

# === Function to check if Mikado preparation exists ===
check_mikado_prepare() {
    if [ -d "working_dir/mikado_op" ]; then
        print_green "Mikado output directory found. Skipping configure and prepare."
        return 0
    else
        print_red "Mikado output directory not found. Running configure and prepare..."
        return 1
    fi
}


# === Function to check if Diamond database exists ===
check_diamond_db() {
    if [ -f "working_dir/prot_evidence.dmnd" ]; then
        print_green "Diamond database found. Skipping creation."
        return 0
    else
        print_red "Diamond database not found. Creating..."
        return 1
    fi
}

# === Function to check if Diamond blastx results exist ===
check_diamond_blastx() {
    if [ -f "working_dir/mikado_prepared.diamond.tsv" ]; then
        print_green "Diamond blastx output found. Skipping blastx."
        return 0
    else
        print_red "Diamond blastx output not found. Running blastx..."
        return 1
    fi
}

# === Function to check if TransDecoder output exists ===
check_transdecoder() {
    if [ -f "working_dir/TransDecoder_op/mikado_prepared.fasta.transdecoder.bed" ]; then
        print_green "TransDecoder output found. Skipping ORF prediction."
        return 0
    else
        print_red "TransDecoder output not found. Running TransDecoder..."
        return 1
    fi
}

# === Function to check if CPC2 result exists ===
check_cpc2() {
    if [ -f "CPC2_result.txt" ] || [ -f "tmp_dir/CPC2_result.txt" ]; then
        print_green "CPC2 result found. Skipping coding potential calculation."
        return 0
    else
        print_red "CPC2 result not found. Running CPC2..."
        return 1
    fi
}

# === Function to check if Salmon quantification exists ===
check_salmon_quant() {
    if [ -f "working_dir/salmon_quant_mikado_prepared/quant.sf" ]; then
        print_green "Salmon quantification found. Skipping quantification."
        return 0
    else
        print_red "Salmon quantification not found. Running Salmon..."
        return 1
    fi
}

# === Function to check additional TPM tmp file ===
check_additional_tmp_tpm() {
    if [[ -f "working_dir/tmp_additional_tpm.tsv" || -f tmp_dir/tmp_additional_tpm.tsv ]]; then
        print_green "Additional temporary TPM file found. Skipping creation."
        return 0
    else
        print_red "Additional temporary TPM file not found. Creating..."
        return 1
    fi
}

# === Function to check additional CPC2 file ===
check_additional_cpc2() {
    if [ -f "working_dir/additional_cpc2.tsv" ]; then
        print_green "Additional CPC2 file found. Skipping creation."
        return 0
    else
        print_red "Additional CPC2 file not found. Creating..."
        return 1
    fi
}

# === Function to check additional TPM file ===
check_additional_tpm() {
    if [ -f "working_dir/additional_tpm.tsv" ]; then
        print_green "Final additional TPM file found. Skipping creation."
        return 0
    else
        print_red "Final additional TPM file not found. Adding missing entries..."
        return 1
    fi
}

# === Function to check external scores file ===
check_external_scores() {
    if [ -f "working_dir/external_scores_tpm_cpc2.tsv" ]; then
        print_green "External scores file found. Skipping creation."
        return 0
    else
        print_red "External scores file not found. Creating..."
        return 1
    fi
}

# === Function to check if Mikado serialise and pick output exists ===
check_mikado_serialise_pick() {
    if [ -f "working_dir/mikado_op/mikado_pick.gff" ]; then
        print_green "Mikado serialise and pick output found. Skipping operation."
        return 0
    else
        print_red "Mikado serialise and pick output not found. Running Mikado..."
        return 1
    fi
}

# === Function to check coding models GFF ===
check_coding_models_gff() {
    if [ -f "working_dir/coding_mikado_models.gff" ]; then
        print_green "Coding models GFF found. Skipping filtering."
        return 0
    else
        print_red "Coding models GFF not found. Filtering..."
        return 1
    fi
}

# === Function to check renamed coding models GFF ===
check_renamed_coding_models() {
    if [ -f "working_dir/renamed_coding_models_mikado.gff" ]; then
        print_green "Renamed coding models GFF found. Skipping renaming."
        return 0
    else
        print_red "Renamed coding models GFF not found. Renaming..."
        return 1
    fi
}

# === Function to check noncoding models GFF ===
check_noncoding_models_gff() {
    if [ -f "working_dir/noncoding_mikado_models.gff" ]; then
        print_green "Noncoding models GFF found. Skipping creation."
        return 0
    else
        print_red "Noncoding models GFF not found. Generating..."
        return 1
    fi
}

# === Function to check extracted Mikado transcripts ===
check_transcripts_mikado_pick() {
    if [ -f "working_dir/transcripts_mikado_pick.fasta" ]; then
        print_green "Mikado transcripts FASTA found. Skipping extraction."
        return 0
    else
        print_red "Mikado transcripts FASTA not found. Extracting..."
        return 1
    fi
}

# === Function to check if GMAP alignment GFF for mikado_pick exists ===
check_gmap_alignment_mikado_pick() {
    if [ -f "working_dir/transcripts_mikado_pick.gff" ]; then
        print_green "GMAP alignment GFF for mikado_pick found. Skipping alignment."
        return 0
    else
        print_red "GMAP alignment GFF for mikado_pick not found. Running alignment..."
        return 1
    fi
}


# === Function to check Mikado GTF ===
check_gmap_alignment_conv_to_gtf() {
    if [ -f "working_dir/transcripts_mikado_pick.gtf" ]; then
        print_green "Mikado transcripts GTF found. Skipping conversion."
        return 0
    else
        print_red "Mikado transcripts GTF not found. Converting..."
        return 1
    fi
}

# === Function to check merged StringTie GTF ===
check_reformat_stringtie_gtf() {
    if [ -f "working_dir/reformat_transcripts_mikado_pick.gtf" ]; then
        print_green "Reformatted Mikado transcripts GTF found. Skipping generation."
        return 0
    else
        print_red "Reformatted Mikado transcripts GTF not found. Reformatting..."
        return 1
    fi
}

# === Function to check Funannotate prediction output ===
check_funannotate_prediction() {
    if [ -f "working_dir/predict_results/predict_misc/evm.round1.gff3.sorted.gff" ]; then
        print_green "Funannotate prediction results found. Skipping prediction."
        return 0
    else
        print_red "Funannotate prediction results not found. Running..."
        return 1
    fi
}

# === Function to check Barrnap GFF ===
check_barrnap() {
    if [ -f "working_dir/reformat_barrnap.gff" ]; then
        print_green "Barrnap GFF found. Skipping annotation."
        return 0
    else
        print_red "Barrnap GFF not found. Running Barrnap..."
        return 1
    fi
}

# === Function to check if EVM output with appended genes exists ===
check_evm_appended() {
    if [ -f "tmp_dir/evm.round1.gff3.sorted.gff" ]; then
        print_green "EVM GFF with appended genes found. Skipping appending."
        return 0
    else
        print_red "EVM GFF with appended genes not found. Appending..."
        return 1
    fi
}

# === Function to check temporary cleaned and renamed GFF ===
check_final_gff() {
    if [ -f "longest_iso_tmp_structural.gff" ] || [ -f "tmp_dir/longest_iso_tmp_structural.gff" ]; then
        print_green "First structural annotation GFF found. Skipping creation."
        return 0
    else
        print_red "First structural annotation GFF not found. Creating..."
        return 1
    fi
}

# === Function to check for discrepancy-free annotation ===
check_discrepancy() {
    if [ -f final_struct_"$species"_"$strain".gff ]; then
        print_green "Fixed "$species" annotation found. Skipping creation."
        return 0
    else
        print_red "Fixed "$species" annotation not found. Creating..."
        return 1
    fi    
}

# === Function to check table2asn has been run ===
check_table2asn() {
    if [ -d table2asn ] || [ -d working_dir/table2asn ]; then
        print_green "Table2asn folder found. Skipping creation."
        return 0
    else
        print_red "Table2asn folder not found. Creating..."
        return 1
    fi    
}

# === Function to check final protein FASTA has been extracted ===
check_final_protein_fasta() {
    if [ -f "final_struct_"$species"_"$strain"_prot.faa" ]; then
        print_green "Final protein FASTA found. Skipping extraction."
        return 0
    else
        print_red "Final protein FASTA not found. Extracting..."
        return 1
    fi
}

# === Function to check ncRNA transcripts FASTA has been extracted ===
check_nc_transcripts_fasta() {
    if [ -f "ncRNA_transcripts_"$species"_"$strain".fasta" ]; then
        print_green "Final ncRNA transcripts FASTA found. Skipping extraction."
        return 0
    else
        print_red "Final ncRNA transcripts FASTA not found. Extracting..."
        return 1
    fi
}

# === Function to check final BUSCO run protein has been run ===
check_busco_prot() {
    if [ -f "busco_prot/short_summary.specific."$busco_lineage"_odb12.busco_prot.txt" ] || [ -f "working_dir/busco_prot/short_summary.specific."$busco_lineage"_odb12.busco_prot.txt" ]; then
        print_green "Final BUSCO run on protein found. Skipping."
        return 0
    else
        print_red "Final BUSCO run on protein not found. Running..."
        return 1
    fi
}

# === Function to check if functional annotation needs to be run or not ===
if ! $no_functional_anno; then
    echo "Functional annotation will be run"
else
    echo "Skipping functional annotation as requested."
fi

# === Function to check if interproscan has been run ===
check_interpro() {
    if [ -f output/input.fa.xml ]; then
     print_green "InterPro prediction found. Skipping ..."
        return 0
    else
        print_red "InterPro prediction not found. Predicting..."
        return 1
    fi
}

# === Function to check if eggnog has been run ===
check_eggnog() {
    if [ -f output/eggnog_"$species"_"$strain".emapper.annotations ]; then
     print_green "EggNOG prediction found. Skipping ..."
        return 0
    else
        print_red "EggNOG prediction not found. Predicting..."
        return 1
    fi
}

# === Function to check if KOfamscan has been run ===
check_kofamscan() {
    if [ -f tmp_dir/filtered_KOfamscan ]; then
     print_green "KOfamscan prediction found. Skipping ..."
        return 0
    else
        print_red "KOfamscan prediction not found. Predicting..."
        return 1
    fi
}

# === Function to check if Infernal has been run ===
check_infernal() {
    if [ -f tmp_dir/Infernal_table.txt ]; then
     print_green "Infernal prediction found. Skipping ..."
        return 0
    else
        print_red "Infernal prediction not found. Predicting..."
        return 1
    fi
}

# === Function to check if Signalp6 has been run ===
check_signalp() {
    if [ -f output/signalp_op_"$species"_"$strain"/prediction_results.txt  ]; then
     print_green "Signalp6 prediction found. Skipping ..."
        return 0
    else
        print_red "Signalp6 prediction not found. Predicting..."
        return 1
    fi
}

# === Function to check if phobius has been run ===
check_phobius() {
    if [ -f output/phobius_"$species"_"$strain".txt ]; then
     print_green "Phobius prediction found. Skipping ..."
        return 0
    else
        print_red "Phobius prediction not found. Predicting..."
        return 1
    fi
}

# === Function to check if tmbed has been run ===
check_tmbed() {
    if [ -f output/Tmbed_"$species"_"$strain".txt ]; then
     print_green "Tmbed prediction found. Skipping ..."
        return 0
    else
        print_red "Tmbed prediction not found. Predicting..."
        return 1
    fi
}


# === Function to check if effectorP3 has been run ===
check_effectorP3() {
    if [ -f output/effectorP3_"$species"_"$strain".txt ]; then
     print_green "EffectorP3 prediction found. Skipping ..."
        return 0
    else
        print_red "EffectorP3 prediction not found. Predicting..."
        return 1
    fi
}


# === Function to check if AntiSMASH has been run ===
check_antismash() {
    if [ -f output/antismash/antismash_"$species"_"$strain".gbk ]; then
     print_green "AntiSMASH prediction found. Skipping ..."
        return 0
    else
        print_red "AntiSMASH prediction not found. Predicting..."
        return 1
    fi
}


# === Function to check if functional_annotation_"$species"_"$strain" has been run ===
check_functional_annotation() {
    if [ -d functional_annotation_"$species"_"$strain" ]; then
     print_green "Functional_annotation folder found. Skipping ..."
        return 0
    else
        print_red "Functional_annotation not found. Appending annotations..."
        return 1
    fi
}

###PIPELINE START
# === Create Folders and agat config===
if ! check_folders_exist; then
    mkdir -p tmp_dir working_dir logs
    agat config --expose
    sed 's/verbose: 1/verbose: 0/g' agat_config.yaml | sed 's/progress_bar: true/progress_bar: false/g' | sed 's/log: true/log: false/g' > tmp && mv tmp agat_config.yaml
fi

# === BUSCO ===
if ! check_busco_results; then
    docker run --rm -v $PWD:/busco_wd/ -v $(realpath "$input_fasta"):/busco_wd/genome.fa ezlabgva/busco:v6.0.0_cv1 bash -c " busco -m genome -f --miniprot -i genome.fa -l "$busco_lineage"_odb12 -o working_dir/busco_genome -c "$threads" " > logs/busco_genome.log
fi

# === Reformat miniprot GFF ===
if ! check_miniprot_results; then
   awk 'BEGIN{FS=OFS="\t"}{if($0!~"#") print}' working_dir/busco_genome/run_"$busco_lineage"*/miniprot_output/*"$busco_lineage"*.gff | sed 's/ /_/g' > tmp_dir/busco_genome.gff 
fi

# === Fix miniprot GFF CDS phases ===
if ! check_miniprot_fix; then
    agat_sp_fix_cds_phases.pl --gff tmp_dir/busco_genome.gff -f $input_fasta -o working_dir/busco.gff 2> logs/agat_sp_fix_cds_phases_busco_genome.log 
fi

# Set flag based on whether Helixer_lineage path is provided
use_helixer=false
if [[ -n "$Helixer_lineage" ]]; then
    use_helixer=true
fi

# Helixer AB Initio prediction
if [[ "$use_helixer" == true ]]; then
    docker run --gpus all --rm \
        -v "$(realpath "$input_fasta")":/home/helixer_user/wd/genome.fa \
        -v "$PWD/":/home/helixer_user/wd \
        gglyptodon/helixer-docker:helixer_v0.3.4_cuda_12.2.2-cudnn8 \
        bash -c "Helixer/scripts/fetch_helixer_models.py && Helixer.py --lineage "$Helixer_lineage" \
        --fasta-path wd/genome.fa --species "$species" --gff-output-path /home/helixer_user/wd/working_dir/Helixer_"$species".gff > logs/helixer_run.log "

else
    echo "--Helixer_lineage flag was not provided, will use the input from --Helixer_gff"
fi

# Set flag based on whether Helixer_lineage path is provided
use_helixer_gff=false
if [[ -n "$Helixer_gff" ]]; then
    use_helixer_gff=true
fi

# === Helixer GFF to GTF conversion ===
if ! check_helixer_gtf; then
    if [[ "$use_helixer_gff" == true ]]; then
        cp "$Helixer_gff" working_dir/Helixer_"$species".gff
        agat_convert_sp_gff2gtf.pl --gff "$Helixer_gff" -o working_dir/helixer.gtf > logs/agat_convert_helixer_gff2gtf.log
    else
        agat_convert_sp_gff2gtf.pl --gff working_dir/Helixer_"$species".gff -o working_dir/helixer.gtf > logs/agat_convert_helixer_gff2gtf.log
    fi
fi

# === STAR genome index ===
if ! check_star_genome_index; then
    STAR --runThreadN "$threads" --runMode genomeGenerate --genomeDir tmp_dir/STAR_index_files --genomeFastaFiles "$input_fasta" --sjdbGTFfile working_dir/helixer.gtf --genomeSAindexNbases 11 > logs/star_index.log
fi

# === STAR mapping ===
if ! check_star_mapping; then
    for i in "${STAR_manifest[@]}"; do
        name=$(basename "$i")
        echo "Mapping reads from $i"
        STAR --runThreadN "$threads" --twopassMode Basic --genomeDir tmp_dir/STAR_index_files --genomeSAindexNbases 11 --readFilesManifest $i --readFilesCommand zcat \
        --outFileNamePrefix tmp_dir/"${name}"_ --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMstrandField intronMotif \
        --limitBAMsortRAM "$limitBAMsortRAM" --outBAMcompression 10 > logs/"${name}"_star_mapping.log
    done
fi

# === Portcullis filtering ===
if ! check_portcullis_filtering; then
    for i in "${STAR_manifest[@]}"; do
        name=$(basename "$i")
        echo "Filtering mappings from ${name}_Aligned.sortedByCoord.out.bam"
        portcullis full -t 13 $input_fasta --orientation $orientation --strandedness $strandedness --max_length $max_intron_length --min_cov 5 tmp_dir/"${name}_Aligned.sortedByCoord.out.bam" -b -o working_dir/"Portcullis_${name}" > logs/portcullis_"${name}_"sj.log 2> logs/portcullis_"${name}_"sj_warnings.log
    done
fi

# === Merge Portcullis junctions ===
if ! check_portcullis_merged; then
    cat working_dir/Portcullis*/3-filt/portcullis_filtered.pass.junctions.bed > working_dir/all_splice_junctions.bed
fi

# === Merge BAM files ===
if ! check_bam_merge; then
    samtools merge -o working_dir/merged_mappings.bam working_dir/Portcullis*/portcullis.filtered.bam -@ $threads 2> logs/samtools_merge.log
fi

# === Aletsch transcript assembly ===
if ! check_aletsch_assembly; then
    mkdir -p gtf profile
    ls -1 working_dir/Portcullis*/*.bam | awk '{print $0, $0".bai", "paired_end"}' > tmp_dir/input_bam_list
    aletsch --profile -i tmp_dir/input_bam_list -p profile > logs/preprocess_aletsch.log
    aletsch -i tmp_dir/input_bam_list -o tmp_dir/raw_aletsch.gtf -p profile -d gtf --max_threads $threads --boost_precision --min_splice_bundary_hits 5 --min_transcript_length_base 200 --min_mapping_quality 20 --output_single_exon_transcripts > logs/aletsch.log
    mv gtf profile tmp_dir
fi

# === StringTie transcript assembly ===
if ! check_stringtie_assembly; then
    stringtie working_dir/merged_mappings.bam -G working_dir/helixer.gtf -p "$threads" -o working_dir/stringtie.gtf -v 2> logs/stringtie.log
fi

# === Reformat Aletsch assembly ===
if ! check_aletsch_merged; then
    stringtie --merge -o working_dir/aletsch.gtf -G tmp_dir/raw_aletsch.gtf -l aletsch tmp_dir/raw_aletsch.gtf -v 2> logs/reformat_aletsch.log
fi

# Set flag based on whether lifted annotation is provided
use_jaccard_clip=false
if [[ -n "$jaccard_clip" ]]; then
    use_jaccard_clip=true
fi

# === Trinity genome-guided transcriptome assembly ===
if ! check_trinity_assembly; then
    if [[ "$use_jaccard_clip" == true ]]; then    
        docker run --rm -v $PWD:/usr/local/src/data trinityrnaseq/trinityrnaseq:2.15.2 bash -c \
        "Trinity --max_memory "$RAM_limit_Trinity" --genome_guided_bam data/working_dir/merged_mappings.bam --jaccard_clip --output data/tmp_dir/op_trinity --full_cleanup --genome_guided_max_intron "$max_intron_length" --CPU "$threads" --SS_lib_type "$orientation" --min_contig_length 200 --grid_node_CPU "$threads" --grid_node_max_memory "$RAM_limit_Trinity" " > logs/Trinity_run.log 
    else
        docker run --rm -v $PWD:/usr/local/src/data trinityrnaseq/trinityrnaseq:2.15.2 bash -c \
        "Trinity --max_memory "$RAM_limit_Trinity" --genome_guided_bam data/working_dir/merged_mappings.bam --output data/tmp_dir/op_trinity --full_cleanup --genome_guided_max_intron "$max_intron_length" --CPU "$threads" --SS_lib_type "$orientation" --min_contig_length 200 --grid_node_CPU "$threads" --grid_node_max_memory "$RAM_limit_Trinity" " > logs/Trinity_run.log
    fi
fi

# === Clean poly-A tails from Trinity transcripts ===
if ! check_seqcleaned_transcripts; then
    seqclean tmp_dir/op_trinity.Trinity-GG.fasta 2> logs/seqclean.log
fi

# === Build GMAP genome index ===
if ! check_gmap_index; then
    gmap_build -d gmap_genome_index -D tmp_dir/ -k 13 "$input_fasta" > logs/gmap_build.log 2> logs/gmap_build_warn.log
fi

# === Align transcripts to genome with GMAP ===
if ! check_gmap_alignment; then
    gmap.sse42 -t "$threads" -n 0 -D ./ -d tmp_dir/gmap_genome_index -f gff3_gene op_trinity.Trinity-GG.fasta.clean > working_dir/trinity.gff 2> logs/gmap_trinity_alignment.log
    mv op_trinity.Trinity-GG.fasta.c* err_seqcl_op_trinity.Trinity-GG.fasta.log outparts_cln.* seqcl_op_trinity.Trinity-GG.fasta.log cleaning_1/ tmp_dir/
fi

# Set flag based on whether lifted annotation is provided
use_lifted_annotation=false
if [[ -n "$lifted_annotation" ]]; then
    use_lifted_annotation=true
fi

# === Reformat lifted annotation for prediction input ===
if [[ "$use_lifted_annotation" == true ]]; then
    if ! check_reformat_lifted_anno; then
        agat_sp_fix_cds_phases.pl --gff "$lifted_annotation" -f "$input_fasta" -o tmp_dir/fix_cds_lifted_anno.gff > logs/fix_cds_lifted_anno.log 2> logs/fix_cds_lifted_anno_warn.log
        awk 'BEGIN{FS=OFS="\t"}{if($3~"mRNA") print $9}' tmp_dir/fix_cds_lifted_anno.gff | sed 's/;/\t/g' | sed 's/Parent=//g' \
        | cut -f2 | grep -w -F -f - tmp_dir/fix_cds_lifted_anno.gff | sed 's/,""//g' | awk 'BEGIN{print "##gff-version 3";FS=OFS="\t"} {print}' > working_dir/liftover.gff
    fi
fi

# === Mikado configure and prepare workflow ===
if ! check_mikado_prepare; then
    cat "$prot_evidence" working_dir/busco_genome/tmp/refseq_db.faa > prot_evidence.faa
    docker run --rm -v $(realpath $mikado_config):/global/mk_config.tsv \
        -v $(realpath $input_fasta):/global/genome.fa \
        -v prot_evidence.faa:/global/prot_evidence.fasta \
        -v $PWD/working_dir:/global/ \
        baderlab/mikado:ubuntu22_mikado2.3.2 bash -c \
        "cp /usr/local/lib/python3.10/dist-packages/Mikado/configuration/scoring_files/HISTORIC/*.yaml /usr/local/lib/python3.10/dist-packages/Mikado/configuration/scoring_files/ ; \
         mikado configure --list mk_config.tsv --max-intron-length "$max_intron_length" --reference genome.fa --strand-specific --mode permissive --check-references \
         --scoring "$scoring" --out-dir mikado_op --junctions all_splice_junctions.bed --blast_targets prot_evidence.fasta --use-transdecoder --threads "$threads" --yaml configuration.yaml; \
         mikado prepare --proc "$threads" --json-conf configuration.yaml"
fi

# === Create Diamond database for protein evidence ===
if ! check_diamond_db; then
    diamond makedb --in prot_evidence.faa -d working_dir/prot_evidence.dmnd 2> logs/diamond_db.log
    rm prot_evidence.faa
fi

# === Run Diamond blastx against protein evidence ===
if ! check_diamond_blastx; then
    diamond blastx --query working_dir/mikado_op/mikado_prepared.fasta \
    --db working_dir/prot_evidence.dmnd --out working_dir/mikado_prepared.diamond.tsv \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop \
    --max-target-seqs 5 --threads $threads 2> logs/diamond_job.log
fi

# === Run TransDecoder for ORF prediction ===
if ! check_transdecoder; then
    TransDecoder.LongOrfs -t working_dir/mikado_op/mikado_prepared.fasta --output_dir working_dir/TransDecoder_op > logs/TransDecoder_longorf.log 2> logs/TransDecoder_longorf_warning.log
    TransDecoder.Predict -t working_dir/mikado_op/mikado_prepared.fasta --output_dir working_dir/TransDecoder_op > logs/TransDecoder_predict.log 2> logs/TransDecoder_predict_warning.log
fi

# === Run CPC2 coding potential calculation ===
if ! check_cpc2; then
    docker run --rm -v $PWD:/data -v $(realpath working_dir/mikado_op/mikado_prepared.fasta):/data/tmp_dir/transcripts.fasta ftricomi/cpc2:latest /CPC2_standalone-1.0.1/bin/CPC2.py -i tmp_dir/transcripts.fasta -o tmp_dir/CPC2_result 2> logs/CPC2.log
fi

# === Quantify transcript abundance with Salmon ===
if ! check_salmon_quant; then
    for i in "${STAR_manifest[@]}"; do cut -f1 "$i" | xargs cat >> left.fq.gz ; done
    for i in "${STAR_manifest[@]}"; do cut -f2 "$i" | xargs cat >> right.fq.gz ; done
     docker run -it --rm -v $PWD:/data quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0 bash -c "salmon index --index data/working_dir/mikado_prepared_index -t data/working_dir/mikado_op/mikado_prepared.fasta -p "$threads" " > logs/salmon_index.log
     docker run -it --rm -v $PWD:/data quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0 bash -c "salmon quant --seqBias --posBias --gcBias -g data/working_dir/mikado_op/mikado_prepared.gtf -p "$threads" -l A  \
     -i data/working_dir/mikado_prepared_index -1 data/left.fq.gz -2 data/right.fq.gz --output data/working_dir/salmon_quant_mikado_prepared " > logs/salmon_quant.log 
     rm left.fq.gz right.fq.gz 
fi

# === Create additional TPM input for Mikado serialise ===
if ! check_additional_tmp_tpm; then
    tail -n +2 working_dir/salmon_quant_mikado_prepared/quant.sf | awk 'BEGIN{FS=OFS="\t"}{print $1,$4}' | sort > working_dir/tmp_additional_tpm.tsv; cp working_dir/tmp_additional_tpm.tsv tmp_dir/tmp_additional_tpm.tsv
fi

# === Create additional CPC2 input for Mikado serialise ===
if ! check_additional_cpc2; then
    tail -n +2 tmp_dir/CPC2_result.txt | awk 'BEGIN{FS=OFS="\t"}{print $1,$7}' | sort > working_dir/additional_cpc2.tsv
fi

# === Add transcripts with zero TPM to additional TPM file ===
if ! check_additional_tpm; then
    cut -f1 working_dir/tmp_additional_tpm.tsv | grep -v -F -w -f - working_dir/additional_cpc2.tsv | awk 'BEGIN{FS=OFS="\t"}{print $1,"0.001"}' | sort >> working_dir/tmp_additional_tpm.tsv; mv working_dir/tmp_additional_tpm.tsv working_dir/additional_tpm.tsv
fi

# === Create combined CPC2 + TPM external scores file ===
if ! check_external_scores; then
    paste working_dir/additional_cpc2.tsv working_dir/additional_tpm.tsv | \
    bedtools groupby -g 1 -c 2,4 -o collapse | \
    sed 's/0\.000000/0.000001/g' | awk 'BEGIN{print "tid\tCPC\ttpm"; FS=OFS="\t"}{print}' > working_dir/external_scores_tpm_cpc2.tsv
fi

# === Mikado serialise and best transcripts picking ===
if ! check_mikado_serialise_pick; then
    docker run --rm -v $(realpath $prot_evidence):/global/prot_evidence.fasta \
        -v $(realpath $input_fasta):/global/genome.fa \
        -v $PWD/working_dir:/global/ \
        baderlab/mikado:ubuntu22_mikado2.3.2 bash -c "cp /usr/local/lib/python3.10/dist-packages/Mikado/configuration/scoring_files/HISTORIC/*.yaml /usr/local/lib/python3.10/dist-packages/Mikado/configuration/scoring_files/ ; mikado serialise --procs "$threads" --json-conf configuration.yaml --orfs TransDecoder_op/mikado_prepared.fasta.transdecoder.bed --transcripts mikado_op/mikado_prepared.fasta --external-scores external_scores_tpm_cpc2.tsv --tsv mikado_prepared.diamond.tsv --blast-targets prot_evidence.fasta --genome genome.fa ; mikado pick --procs "$threads" --json-conf configuration.yaml --loci-out mikado_pick.gff --monoloci-out monoloci_mikado_pick.gff --genome genome.fa"
fi

# === Filter out bad coding gene models and clean GFF for training ===
if ! check_coding_models_gff; then
    awk 'BEGIN{FS=OFS="\t"}{if(($3=="ncRNA") || ($3=="mRNA" && $9~"has_start_codon=False" || $9~"has_stop_codon=False")) print $9}' working_dir/mikado_op/mikado_pick.gff | \
    sed 's/;/\t/g' | sed 's/Parent=//g' | cut -f2 - | grep -w -F -f - working_dir/mikado_op/mikado_pick.gff | \
    awk 'BEGIN{print "##gff-version 3";FS=OFS="\t"}{if($9!~"has_start_codon=True" || $9!~"has_stop_codon=True") print}' | \
    awk 'BEGIN{FS=OFS="\t"}{if($3!~"UTR" && $3!~"CDS") print}' | \
    grep -w -F -v -f - working_dir/mikado_op/mikado_pick.gff | sed 's/Name.*//g' | \
    awk 'BEGIN{FS=OFS="\t"}{if($3!="superlocus" && $3!~"UTR")print}' > working_dir/coding_mikado_models.gff
fi

# === Rename and clean Mikado GFF coding for training ===
if ! check_renamed_coding_models; then
    docker run --rm -v $(realpath $input_fasta):/data/genome.fa \
        -v $PWD/working_dir:/data nextgenusfs/funannotate:v1.8.17 bash -c \
        "funannotate gff-rename -g data/coding_mikado_models.gff -f data/genome.fa -o data/renamed_coding_models_mikado.gff"
fi

# === Write noncoding GFF from total ===
if ! check_noncoding_models_gff; then
    awk 'BEGIN{FS=OFS="\t"}{if(($3=="ncRNA") || ($3=="mRNA" && $9~"has_start_codon=False" || $9~"has_stop_codon=False")) print $9}' working_dir/mikado_op/mikado_pick.gff | \
    sed 's/;/\t/g' | sed 's/Parent=//g' | cut -f2 - | grep -w -F -f - working_dir/mikado_op/mikado_pick.gff | \
    awk 'BEGIN{print "##gff-version 3";FS=OFS="\t"}{if($9!~"has_start_codon=True" || $9!~"has_stop_codon=True") print}' | \
    awk 'BEGIN{FS=OFS="\t"}{if($3!~"UTR" && $3!~"CDS") print}'  | sed 's/Name=.*//g' | sed 's/mRNA/ncRNA/g' | \
    awk 'BEGIN{FS=OFS="\t"}{if($3!~"UTR" && $3!~"CDS") print}' | sed 's/alias.*//g' | sed 's/ncRNA_gene/gene/g' | sed 's/fpkm=.*//g' > working_dir/noncoding_mikado_models.gff
fi

# === Extract transcripts from Mikado overall transcript assembly-annotation ===
if ! check_transcripts_mikado_pick; then
    agat_sp_extract_sequences.pl --type exon --merge -g working_dir/mikado_op/mikado_pick.gff -f $input_fasta -o working_dir/transcripts_mikado_pick.fasta > logs/agat_extract_seq.log
fi

# === Align mikado_pick_transcripts to genome and create gff ===
if ! check_gmap_alignment_mikado_pick; then
    gmap.sse42 -t $threads -n 0 -D ./ -d tmp_dir/gmap_genome_index -f gff3_gene working_dir/transcripts_mikado_pick.fasta > working_dir/transcripts_mikado_pick.gff 2> logs/transcripts_mikado_pick_gmap_alignment.log
fi

# === Convert Mikado GFF to GTF ===
if ! check_gmap_alignment_conv_to_gtf; then
    agat_convert_sp_gff2gtf.pl --gff working_dir/transcripts_mikado_pick.gff -o working_dir/transcripts_mikado_pick.gtf > logs/convert_transcripts_mikado.log
fi

# === Merge StringTie GTF ===
if ! check_reformat_stringtie_gtf; then
    stringtie --merge -G working_dir/transcripts_mikado_pick.gtf -l mikado working_dir/transcripts_mikado_pick.gtf -o working_dir/reformat_transcripts_mikado_pick.gtf -v 2> logs/reformat_mikado_transcripts.log
fi

# Build --other_gff argument based on whether lifted annotation is used
other_gff_arg="--other_gff data/Helixer_"$species".gff:4"
if [[ "$use_lifted_annotation" == true ]]; then
    other_gff_arg="$other_gff_arg data/liftover.gff:2"
fi

# Set flag based on whether GeneMark path is provided
use_genemark=false
if [[ -n "$GeneMark_PATH" ]]; then
    use_genemark=true
fi

# === Run Funannotate ===

if ! check_funannotate_prediction; then
    if [[ "$use_genemark" == true ]]; then
        docker run --rm \
            -v "$(realpath $prot_evidence)":/data/prot_evidence.fasta \
            -v "$(realpath $input_fasta)":/data/genome.fa \
            -v "$PWD/working_dir":/data \
            -v "$GeneMark_PATH"/:/GeneMark-ETP \
            nextgenusfs/funannotate:v1.8.17 bash -c \
            " funannotate predict -i data/genome.fa -o data/predict_results -s \"$species\" \
            --GENEMARK_PATH /GeneMark-ETP/bin/gmes/ --genemark_mode ET \
            --pasa_gff data/renamed_coding_models_mikado.gff \
            --protein_evidence data/prot_evidence.fasta \
            $other_gff_arg \
            --rna_bam data/merged_mappings.bam \
            --stringtie data/reformat_transcripts_mikado_pick.gtf \
            --cpus $threads --optimize_augustus \
            --weights augustus:1 pasa:5 hiq:3 codingquarry:1 glimmerhmm:1 snap:1 genemark:1 \
            --transcript_alignments data/transcripts_mikado_pick.gff"
    else
        docker run --rm \
            -v "$(realpath $prot_evidence)":/data/prot_evidence.fasta \
            -v "$(realpath $input_fasta)":/data/genome.fa \
            -v "$PWD/working_dir":/data \
            nextgenusfs/funannotate:v1.8.17 bash -c \
            " funannotate predict -i data/genome.fa -o data/predict_results -s \"$species\" \
            --pasa_gff data/renamed_coding_models_mikado.gff \
            --protein_evidence data/prot_evidence.fasta \
            $other_gff_arg \
            --rna_bam data/merged_mappings.bam \
            --stringtie data/reformat_transcripts_mikado_pick.gtf \
            --cpus $threads --optimize_augustus \
            --weights augustus:1 pasa:5 hiq:3 codingquarry:1 glimmerhmm:1 snap:1 genemark:1 \
            --transcript_alignments data/transcripts_mikado_pick.gff" 
    fi
fi    

# === Run Barrnap ===
if ! check_barrnap; then
    barrnap --kingdom euk --threads $threads "$input_fasta" --reject 0.50 > tmp_dir/raw_barrnap.gff 2> logs/barrnap.log
    awk 'BEGIN{FS=OFS="\t"}{if($9!~"partial") print}' tmp_dir/raw_barrnap.gff > tmp && mv tmp tmp_dir/raw_barrnap.gff
    agat_sp_fix_overlaping_genes.pl -f tmp_dir/raw_barrnap.gff -o working_dir/reformat_barrnap.gff > logs/fix_barrnap.log 2> logs/fix_barrnap_anno_warn.log
fi

# === Append missing gene models ===
if ! check_evm_appended; then
    set +e
    cp working_dir/predict_results/predict_misc/evm.round1.gff3.sorted.gff ./

    for i in working_dir/reformat_barrnap.gff working_dir/predict_results/predict_misc/pasa_predictions.gff3 working_dir/predict_results/predict_misc/other1_predictions.gff3 working_dir/noncoding_mikado_models.gff working_dir/predict_results/predict_misc/genemark.gff working_dir/predict_results/predict_misc/augustus.gff3 working_dir/predict_results/predict_misc/glimmerhmm-predictions.gff3 working_dir/predict_results/predict_misc/snap-predictions.gff3 working_dir/predict_results/predict_misc/trnascan.gff3.sorted.gff3 working_dir/predict_results/predict_misc/other2_predictions.gff3; do
        if [ -f "$i" ]; then
            echo "Integrating gene models from "$i" to final annotation... "
            bedtools intersect -a "$i" -b evm.round1.gff3.sorted.gff -wa -wb -s -v | awk 'BEGIN{FS=OFS="\t"}{if($3=="mRNA" || $3=="tRNA" || $3=="rRNA" || $3=="gene") print $9}' | \
            sed 's/ID=//g' | sed 's/;/\t/g' | cut -f1 | grep -v "^$" | grep -F -w -f - "$i" | sed 's/Alias.*//g' | sed 's/Name.*\t//g' | \
            awk 'BEGIN{FS=OFS}{if($0!~"#" && $3!~"intron" && $3!~"codon" && $3!~"UTR") print}' | sed 's/transcript/mRNA/g' | sed 's/other_pred1/Helixer/g' | sed 's/pasa/Mikado_loci/g' | \
            sed 's/AGAT/Barrnap/g' | sed 's/\[\]/hypothetical protein/g' >> evm.round1.gff3.sorted.gff
        else
            echo "[WARNING] File $i not found, skipping."
        fi
    done
        mv evm.round1.gff3.sorted.gff tmp_dir/
        set -e
fi

# === Clean, rename, and sort GFF ===
if ! check_final_gff; then
    agat_sp_fix_overlaping_genes.pl -g tmp_dir/evm.round1.gff3.sorted.gff -o tmp_dir/fix_olp.gff > logs/fix_olp.log
    agat_sp_fix_features_locations_duplicated.pl -g tmp_dir/fix_olp.gff -o tmp_dir/fix_locdup.gff > logs/fix_locdup.log
    agat_sp_filter_incomplete_gene_coding_models.pl -f "$input_fasta" -g tmp_dir/fix_locdup.gff -o tmp_dir/filter_incomplete_gene.gff > logs/filter_incomplete_gene.gff
    agat_sp_filter_by_ORF_size.pl -s 49 --gff tmp_dir/filter_incomplete_gene.gff -o tmp_dir/orf_filt_incomplete.gff > logs/orf_filt_incomplete.log
    awk '$3=="gene"{match($0,/ID=([^;]+)/,a); print a[1],$5-$4+1}' tmp_dir/orf_filt_incomplete_sup49.gff | sort -k2,2nr | awk -v max="$max_gene_length" '$2 >= max {print}' | cut -f1 -d " " > tmp_dir/kill_list.txt
    agat_sp_filter_feature_from_kill_list.pl --gff tmp_dir/orf_filt_incomplete_sup49.gff --kill_list tmp_dir/kill_list.txt -o tmp_dir/genes_kept1.gff > logs/filter_from_kill_list.log

    docker run --rm \
        -v $(realpath $input_fasta):/data/genome.fa \
        -v $PWD/tmp_dir:/data \
        nextgenusfs/funannotate:v1.8.17 bash -c \
        "funannotate gff-rename -g data/genes_kept1.gff -f data/genome.fa -l $locus_tag -o data/renamed_all.gff"

    awk 'BEGIN{FS=OFS="\t"}{if($3=="ncRNA" && $5-$4+1>=200){gsub(/product=None/, "product=hypothetical lncRNA;ncRNA_class=lncRNA", $0)} else if ($3=="ncRNA" && $5-$4+1<200){gsub(/product=None/, "product=hypothetical ncRNA;ncRNA_class=other", $0)}print}' tmp_dir/renamed_all.gff | \
    sed 's/Alias=.*//g' | sed 's/Name=.*//g' > tmp_dir/tmp_sorted_structural.gff
    agat_sp_keep_longest_isoform.pl --gff tmp_dir/tmp_sorted_structural.gff -o tmp_dir/longest_iso_tmp_structural.gff > logs/longest_iso_tmp_structural.log
fi

# === Remove discrepant gene models and create gbk file ===
if ! check_discrepancy; then
    docker run -v $(realpath $input_fasta):/data/genome.fa -v $PWD/:/data --rm nextgenusfs/funannotate:v1.8.17 bash -c \
    "funannotate gff2tbl -g data/tmp_dir/longest_iso_tmp_structural.gff -f data/genome.fa > data/tmp_dir/longest_iso_tmp_structural.tbl 
    funannotate tbl2gbk -i data/tmp_dir/longest_iso_tmp_structural.tbl -f data/genome.fa -s "$species" -o data/tmp_dir/longest_iso_tmp_structural > data/tmp_dir/genes_to_fix_or_remove 
    tail -n +3 data/tmp_dir/genes_to_fix_or_remove | cut -f1 | grep -v -F -w -f - data/tmp_dir/longest_iso_tmp_structural.gff > data/tmp_dir/genes_kept2.gff
    funannotate gff-rename -g data/tmp_dir/genes_kept2.gff -f data/genome.fa -l $locus_tag -o data/tmp_dir/ren2.gff "
    sed 's/Alias=.*//g' tmp_dir/ren2.gff | awk 'BEGIN{FS=OFS="\t"}{if($3=="ncRNA" && $5-$4+1>=200){gsub(/product=hypothetical lncRNA/, "product=hypothetical lncRNA;ncRNA_class=lncRNA", $0)} else if ($3=="ncRNA" && $5-$4+1<200){gsub(/product=hypothetical ncRNA/, "product=hypothetical ncRNA;ncRNA_class=other", $0)}print}' \
    > final_struct_"$species"_"$strain".gff 
fi

# === Run Table2asn and create gbk file for AntiSMASH ===
if ! check_table2asn; then 
    table2asn -V vb -M n -J -c ewf -euk -t "$submission_template" -gaps-min 10 -l paired-ends -locus-tag-prefix "$locus_tag" \
    -j "[organism="$species"] [strain="$strain"] [gcode="$codon_table"]" -i "$input_fasta" -f final_struct_"$species"_"$strain".gff -Z -outdir working_dir/table2asn -verbose 2> logs/table2asn.log
fi

# === Extract final protein sequences ===
if ! check_final_protein_fasta; then
    agat_sp_extract_sequences.pl --clean_internal_stop --clean_final_stop --codon $codon_table --aa \
    --gff final_struct_"$species"_"$strain".gff -f $input_fasta -o final_struct_"$species"_"$strain"_prot.faa > logs/final_protein_fasta.log
fi

# === Extract noncoding transcripts for infernal prediction ===
if ! check_nc_transcripts_fasta; then
    agat_sp_extract_sequences.pl --type ncRNA --merge -f "$input_fasta" -g final_struct_"$species"_"$strain".gff -o ncRNA_transcripts_"$species"_"$strain".fasta > logs/check_nc_transcripts.log
fi

# === Run BUSCO on proteins ===
if ! check_busco_prot; then
    docker run --rm -v $PWD:/busco_wd/ ezlabgva/busco:v6.0.0_cv1 bash -c " busco -m proteins -i final_struct_"$species"_"$strain"_prot.faa -l "$busco_lineage"_odb12 -o working_dir/busco_prot -f -c "$threads" " > logs/busco_prot.log
fi

# === FUNCTIONAL ANNOTATION ===
if [ "$no_functional_anno" != true ]; then
        mkdir -p output
        echo "Functional annotation will be run in this order : InterProScan - EggNOG - KOfamscan - Infernal - Signalp6 - Phobius - Tmbed - AntiSMASH"
    fi

# === Run InterProScan ===
    if ! check_interpro; then
        docker run --rm -v $PWD/logs/:/logs -v $(realpath $databases/interproscan-5.75-106.0/data/):/opt/interproscan/data -v $PWD/output/:/output -v $PWD/final_struct_"$species"_"$strain"_prot.faa:/input/input.fa interpro/interproscan:5.75-106.0 \
        --input /input/input.fa --formats XML --disable-precalc --output-dir /output/ --cpu 6
    fi

# === Run EggNog ===
    if ! check_eggnog; then
        docker run --rm -v $PWD/logs/:/logs -v $(realpath $databases):/databases -v $PWD/final_struct_"$species"_"$strain"_prot.faa:/media/prot.faa -v $PWD/output/:/output quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_2 bash -c \
        " emapper.py --cpu "$threads" --go_evidence all -m diamond --data_dir databases/ -i media/prot.faa -o output/eggnog_"$species"_"$strain" " > logs/eggnog.log 2> logs/eggnog_warn.log
    fi

# === Run Kofamscan === 
    if ! check_kofamscan; then
        docker run --rm -v $(realpath $databases):/databases -v $PWD/final_struct_"$species"_"$strain"_prot.faa:/media/prot.faa -v $PWD:/data quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2 bash -c " exec_annotation -p databases/profiles/ -k databases/ko_list --cpu "$threads" -f detail-tsv -o data/tmp_dir/raw_KOfamscan_result media/prot.faa "
        bit-filter-KOFamScan-results -i tmp_dir/raw_KOfamscan_result -o tmp_dir/filtered_KOfamscan 
        awk 'BEGIN{FS=OFS="\t"} NR>1 {gsub(/\[EC[:]?/, "\tEC_number", $0); gsub(/\]$/, "", $0); if($2!="NA" && $3!="NA"){print $1,"note","KEGG:"$2; if($4!=""){ec=$4; gsub(/ /,",",ec); print $1,"EC_number",ec} print $1,"product",$3}}' tmp_dir/filtered_KOfamscan | sed 's/EC_number\tEC_number/EC_number\t/g' > working_dir/additional_annos.tsv
    fi

# === Run Infernal ===
    if ! check_infernal; then
        cmscan --cpu "$threads" --cut_ga --rfam --nohmmonly --clanin "$databases"/Rfam.clanin --oskip --fmt 2 --tblout tmp_dir/Infernal_table.txt "$databases"/Rfam.cm ncRNA_transcripts_"$species"_"$strain".fasta > logs/infernal.log
        awk '!/^#|^-/ && NF >= 4 {print $4 "\tdb_xref\t" $3; print $4 "\tncRNA_class\t" $2}' tmp_dir/Infernal_table.txt >> working_dir/additional_annos.tsv
    fi

# === Run SignalP6 ===
    if ! check_signalp; then
        signalp6 --fastafile final_struct_"$species"_"$strain"_prot.faa --organism eukarya --format txt --output_dir output/signalp_op_"$species"_"$strain" --mode slow-sequential -md "$databases" 2> logs/signalp6.log
    fi

# === Run Phobius ===
    if ! check_phobius; then
        "$phobius_PATH"/phobius.pl final_struct_"$species"_"$strain"_prot.faa -short > output/phobius_"$species"_"$strain".txt 2> logs/phobius.log
    fi

# === Run Tmbed ===
    if [ "$no_tmbed" != true ]; then
        if ! check_tmbed; then
          docker run --rm --gpus all -v $(realpath "$databases"):/databases -v $PWD:/data -v $PWD/output/:/output interpro/tmbed:1.0.2 bash -c \
          "tmbed predict --out-format=0 -f data/final_struct_"$species"_"$strain"_prot.faa -m databases/"$databases" -p output/Tmbed_"$species"_$strain".txt > data/logs/Tmbed.log"
          awk 'NR%3==1{id=substr($1,2)} NR%3==0{f=$0;cat="";if(f~/^[.S]+$/&&f~/S/)cat="TMbed:Secreted";else if(f~/[hHbB]/)cat="Tmbed:Transmembrane";if(cat!="")printf "%s\tnote\t%s\n",id,cat}' output/Tmbed_"$species"_$strain".txt >> working_dir/additional_annos.tsv
        fi
    fi

# === Run EffectorP3 ===
if [ "$no_effectorp3" != true ]; then
    if ! check_effectorP3; then
        docker run --rm -v $PWD/output/:/output -v $PWD:/data gabrielerigano/effectorp3:latest bash -c "python EffectorP-3.0/EffectorP.py -f -i /data/final_struct_"$species"_"$strain"_prot.faa -o /output/effectorP3_"$species"_"$strain".txt > /data/logs/Effectorp3.log"
        awk 'BEGIN{FS=OFS="\t"}{if($4=="-") print $1,"note",$5}' output/effectorP3_"$species"_"$strain".txt | tail -n +2 >> working_dir/additional_annos.tsv
    fi
fi

# === Run AntiSMASH === 
    if [ "$no_antismash" != true ]; then
        if ! check_antismash; then
            docker run -v $(realpath working_dir/table2asn/*.gbf):/input/table2asn.gbk -v $PWD/output/:/output --rm antismash/standalone:8.0.4 table2asn.gbk --output-dir /output/antismash \
            --output-basename antismash_"$species"_"$strain" --rre --genefinding-tool none -t fungi --cpus "$threads" --verbose 2> logs/antismash.log
        fi
    fi

# === Run funannotate annotate ===
    if [ "$no_antismash" != true ]; then
        if ! check_functional_annotation; then
            docker run -v $(realpath "$submission_template"):/data/tmp_dir/submission.sbt -v $(realpath "$input_fasta"):/data/genome.fa -v $PWD/:/data --rm nextgenusfs/funannotate:v1.8.17 bash -c \
            "sed s'/1\.4/4\.0/g' venv/lib/python3.8/site-packages/funannotate/downloads.json > tmp/tmp && mv tmp/tmp venv/lib/python3.8/site-packages/funannotate/downloads.json; \
            funannotate setup -i merops uniprot dbCAN pfam go mibig interpro gene2product -u -d /opt/databases/ ; funannotate annotate --gff data/final_struct_"$species"_"$strain".gff \
            --fasta data/genome.fa -o data/functional_annotation_"$species"_"$strain" --cpus "$threads" --annotations data/working_dir/additional_annos.tsv --eggnog data/output/eggnog_"$species"_"$strain".emapper.annotations \
            --iprscan data/output/input.fa.xml --phobius data/output/phobius_"$species"_"$strain".txt --signalp data/output/signalp_op_"$species"_"$strain"/prediction_results.txt --strain "$strain" \
            --sbt data/tmp_dir/submission.sbt --species "$species" " 
        else
            docker run -v $(realpath "$submission_template"):/data/tmp_dir/submission.sbt -v $(realpath "$input_fasta"):/data/genome.fa -v $PWD/:/data --rm nextgenusfs/funannotate:v1.8.17 bash -c \
            "sed s'/1\.4/4\.0/g' venv/lib/python3.8/site-packages/funannotate/downloads.json > tmp/tmp && mv tmp/tmp venv/lib/python3.8/site-packages/funannotate/downloads.json ; \
            funannotate setup -i merops uniprot dbCAN pfam go mibig interpro gene2product -u -d /opt/databases/ ; funannotate annotate --gff data/final_struct_"$species"_"$strain".gff \
            --fasta data/genome.fa -o data/functional_annotation_"$species"_"$strain" --cpus "$threads" --annotations data/working_dir/additional_annos.tsv --eggnog data/output/eggnog_"$species"_"$strain".emapper.annotations \
            --antismash data/output/antismash/antismash_"$species"_"$strain".gbk --iprscan data/output/input.fa.xml --phobius data/output/phobius_"$species"_"$strain".txt \
            --signalp data/output/signalp_op_"$species"_"$strain"/prediction_results.txt --strain "$strain" \
            --sbt data/tmp_dir/submission.sbt --species "$species" "
        fi
            echo -e "\e[32m[SUCCESS]\e[0m bagRNA has completed successfully."
        exit 0
else
        echo -e "Skipping functional annotation (--no_functional_anno flag set), bagRNA has completed successfully"
        exit 0
fi
