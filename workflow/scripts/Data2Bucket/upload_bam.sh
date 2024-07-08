#!/bin/bash
umask 000

# Input parameters
bcl=$1
count_folder=$2
log_folder=$3
pkg_path=$4
bucket_name=$5

# Load packages
gsutil_pkg="$pkg_path/google-cloud-sdk/bin/gsutil"
gcloud_pkg="$pkg_path/google-cloud-sdk/bin/gcloud"

# upload bam files
if [ ! -z "$(ls -A "$count_folder")" ]; then
    failed_uploads="$log_folder/upload/bam_failed_uploads.txt"
    files_to_upload="$log_folder/upload/bam_files_to_upload.txt"
    touch "$failed_uploads" >/dev/null 2>&1
    touch "$files_to_upload" >/dev/null 2>&1
    > "$failed_uploads"

    count_file_list=()
    while IFS= read -r directory; do
        if [[ -f "$directory/outs/possorted_genome_bam.bam" ]]; then
            count_file_list+=("$directory/outs/possorted_genome_bam.bam")
        fi
    done < <(find "$count_folder" -mindepth 1 -maxdepth 1 -type d -print)

    if [ ${#count_file_list[@]} -eq 0 ]; then
        echo "No bam files found."
        exit 0
    fi
    echo "Bam files found:"
    for file in "${count_file_list[@]}"; do
        echo "$file"
    done
    echo "${count_file_list[@]}" > "$files_to_upload"

    for file in "${count_file_list[@]}"; do
        IFS='/' read -ra ADDR <<< "$file"
        sample="${ADDR[-3]}"
        gs_path="gs://$bucket_name/bam/$bcl/$sample/$(basename "$file")"

        if ! "$gsutil_pkg" -q stat "$gs_path"; then
            attempt=0
            success=0
            while [ $attempt -lt 3 ] && [ $success -eq 0 ]; do
                echo "Attempting to upload $file to $gs_path (Attempt $((attempt + 1)))"
                if "$gcloud_pkg" alpha storage cp "$file" "$gs_path"; then
                    success=1
                else
                    echo -e "\nFailed to upload: $file to $gs_path" >&2
                    attempt=$((attempt + 1))
                    sleep 5
                fi
            done
            
            if [ $success -eq 0 ]; then
                echo "$file" >> "$failed_uploads"
            fi
        else
            echo "File already exists and is assumed to be complete."
        fi
    done

    if [ -s "$failed_uploads" ]; then
        echo -e "\nFailed to upload some files, please check $failed_uploads for details." >&2
    fi

    echo "Finished uploading bam files."
else
    echo "No bam folders found."
fi
