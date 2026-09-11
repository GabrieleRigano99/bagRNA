#!/usr/bin/env bash

# === Default Settings ===
database_folder="databases"
force_download=0

# === Color Output Functions ===
print_green()  { echo -e "\e[32m$1\e[0m"; }
print_red()    { echo -e "\e[31m$1\e[0m"; }
print_yellow() { echo -e "\e[33m$1\e[0m"; }

# === Usage Message ===
usage() {
    echo -e "\nUsage: $0 [OPTIONS]"
    echo "Download InterPro, EggNOG, KEGG KO, and RFAM databases into a specified folder."
    echo ""
    echo "Options:"
    echo "  -h, --help                    Show this help message and exit"
    echo "  --database_folder <path>     Path to download the databases (default: ./databases)"
    echo "  -f, --force                   Force redownload even if the folder exists"
    echo ""
    exit 1
}

# === Parse Arguments ===
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        --database_folder)
            if [[ -z "$2" || "$2" == --* ]]; then
                echo "Error: --database_folder requires a path."
                usage
            fi
            database_folder="$2"
            shift 2
            ;;
        -f|--force)
            force_download=1
            shift
            ;;
        *)
            echo "Error: Unknown argument '$1'"
            usage
            ;;
    esac
done

# === Check if database folder exists ===
check_dbs_exist() {
    [ -d "$database_folder" ]
}

# === Validate Presence of Expected Files ===
validate_dbs() {
    local missing=0

    [[ -d "$database_folder/eggnog_proteins" ]] || { print_yellow "⚠️ Missing: EggNOG data"; missing=1; }
    [[ -d "$database_folder/interproscan" || -f "$database_folder/interproscan.xml" ]] || { print_yellow "⚠️ Missing: InterProScan data"; missing=1; }
    [[ -f "$database_folder/Rfam.cm" ]] || { print_yellow "⚠️ Missing: Rfam.cm"; missing=1; }
    [[ -f "$database_folder/Rfam.clanin" ]] || { print_yellow "⚠️ Missing: Rfam.clanin"; missing=1; }
    [[ -d "$database_folder/profiles" ]] || { print_yellow "⚠️ Missing: KEGG KO profiles"; missing=1; }
    [[ -f "$database_folder/ko_list" ]] || { print_yellow "⚠️ Missing: KEGG ko_list"; missing=1; }

    return $missing
}

# === Download All Databases ===
download_databases() {
    mkdir -p "$database_folder"

    print_green "⬇️ Downloading EggNOG..."
    download_eggnog_data.py -F -P --data_dir "$database_folder/"

    print_green "⬇️ Downloading InterProScan..."
    curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/alt/interproscan-data-5.75-106.0.tar.gz \
        && tar -xzf interproscan-data-5.75-106.0.tar.gz -C "$database_folder" \
        && rm interproscan-data-5.75-106.0.tar.gz

    print_green "⬇️ Downloading Rfam..."
    wget -q -P "$database_folder" https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz \
        && gzip -df "$database_folder/Rfam.cm.gz"
    wget -q -P "$database_folder" https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
    cmpress "$database_folder"/Rfam.cm

    print_green "⬇️ Downloading KEGG KO..."
    curl -s -L -O ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz \
        && tar -xzf profiles.tar.gz -C "$database_folder" \
        && rm profiles.tar.gz
    curl -s -L -o "$database_folder/ko_list.gz" ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz \
        && gzip -df "$database_folder/ko_list.gz"
}

# === Main Execution ===
if check_dbs_exist && [[ "$force_download" -eq 0 ]]; then
    print_green "✅ Database folder already exists: $database_folder"
    print_green "🔍 Validating contents..."
    if validate_dbs; then
        print_green "✅ All required database files are present and correct."
        exit 0
    else
        print_red "⚠️ Some database files are missing or incomplete."
        print_red "💡 Re-run with --force to redownload all."
        exit 1
    fi
fi

if [[ "$force_download" -eq 1 ]]; then
    print_yellow "⚠️ --force specified. Redownloading all databases..."
fi

download_databases
print_green "✅ All databases have been downloaded into: $database_folder"
