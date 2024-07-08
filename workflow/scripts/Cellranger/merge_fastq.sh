#!/bin/bash
umask 000

base_sheet_folder=$1
link_file_type=$2
download_files=$3
use_cluster=$4

current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source "$current_dir/../../config/config.sh"
merge_csv="$base_sheet_folder/merge.csv"
base_path="${BASE_DATA_PATH}"
PKG_PATH="${PKG_PATH}"
GOOGLE_CLOUD_BUCKET="${GOOGLE_CLOUD_BUCKET}"

# add index to avoid overwriting for the same sample
increment_suffix() {
    filename=$1
    base_name=$(echo "$filename" | grep -oP '.*(?=_\d{3}\.fastq\.gz)')
    suffix=$(echo "$filename" | grep -oP '(?<=_)\d{3}(?=\.fastq\.gz)')
    new_suffix=$(printf "%03d" $((10#$suffix + 1)))
    echo "${base_name}_${new_suffix}.fastq.gz"
}

# cat fastq files with the same file type
merge_files() {
    file_type=$1
    folder_path="$to_path/before"
    files=($(ls ${folder_path}/*_S*_L*_${file_type}_*.fastq.gz 2> /dev/null | sort))
    if [ ${#files[@]} -gt 0 ]; then
        cat "${files[@]}" > "${to_path}/$(basename ${files[0]})"
    fi
}
# check and merge fastq files
check_and_merge() {
    file_type=$1
    if ! ls ${to_path}/*_S*_L*_${file_type}_*.fastq.gz 2> /dev/null; then
        merge_files "$file_type"
    fi
}

# link fastqs from the additional BCL
link_fastqs() {
    bcl=$1
    add_index=$2
    index_nd=$3

    if [ ! -z "$add_index" ]; then
        if [ "$add_index" == "$index_nd" ]; then
            index="$index_nd"
        else
            index="$add_index"
        fi
        from_path="$base_path/$bcl/fastq/$index"
    else
        index="$index_nd"
        from_path="$base_path/$bcl/fastq/$index"
    fi

    # download files from google cloud if files do not exist or forced to download
    file_exist=false
    if ls ${from_path}/${index}*_S*_L*_R1* 1> /dev/null 2>&1 || ls ${from_path}/${index}*_S*_L*_R2* 1> /dev/null 2>&1; then
        file_exist=true
    fi
    if [ "$download_files" == "true" ] || [ "$file_exist" == "false" ]; then
        base_fastq_path="$base_path/$bcl/fastq"
        base_log_path="$base_path/$bcl/log"
        download_files "$bcl" "$base_fastq_path" "" "$index" "$base_log_path" "$PKG_PATH" "$GOOGLE_CLOUD_BUCKET" "$use_cluster"
    fi

    to_path="$base_path/$main_bcl/fastq/$index_nd"
    mkdir -p "$to_path" 2> /dev/null

    found=0
    for file in ${from_path}/*_S*_L*_*.fastq.gz; do
        if [[ -f "$file" ]]; then
            found=1
            base_name=$(basename "$file")
            if [[ "$add_index" != "$index_nd" ]]; then
                base_name="${base_name/$add_index/$index_nd}"
            fi
            sample_prefix=$(echo $base_name | grep -oP '^[^_]+')
            new_filename=$(increment_suffix "$base_name" | sed "s/^${sample_prefix}_//")
            ln -sf "$file" "${to_path}/${sample_prefix}_${new_filename}"
        fi
    done
    if [ $found -eq 0 ]; then
        echo "Error: No matching files found in $from_path."
        exit 1
    fi
}

if [ -f "$merge_csv" ]; then
    echo -e "\nLink fastqs from additional BCLs..."
    tail -n +2 "$merge_csv" | while IFS=',' read -r main_bcl rna_bcl add_rna_index sb_bcl add_sb_index add_sb_puck name rna_index_nd sb_index_nd; do
        
        # Merge fastq files from addtional RNA BCL
        if [ "$link_file_type" == "rna" ]; then
            if [ ! -z "$rna_bcl" ]; then
                echo "Link fastqs for RNA Index: $rna_index_nd from $rna_bcl to $main_bcl"
                link_fastqs "$rna_bcl" "$add_rna_index" "$rna_index_nd"
            fi
        fi

        # Merge fastq files from addtional SB BCL
        if [ "$link_file_type" == "sb" ]; then
            if [ ! -z "$sb_bcl" ]; then
                echo "Link fastqs for SB Index: $sb_index_nd from $sb_bcl to $main_bcl"
                link_fastqs "$sb_bcl" "$add_sb_index" "$sb_index_nd"
            fi
        fi

    done
fi