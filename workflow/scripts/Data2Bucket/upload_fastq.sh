#!/bin/bash
umask 000

# Input parameters
bcl=$1
fastq_folder=$2
log_folder=$3
pkg_path=$4
bucket_name=$5

# Load packages
gsutil_pkg="$pkg_path/google-cloud-sdk/bin/gsutil"
gcloud_pkg="$pkg_path/google-cloud-sdk/bin/gcloud"

# upload fastq files
if [ ! -z "$(ls -A "$fastq_folder")" ]; then
    fastq_folder_list=$(find "$fastq_folder" -mindepth 1 -maxdepth 1 -type d ! -lname '*' ! -name 'mkfastq*' -print)
    echo -e "Fastq folders found: \n$fastq_folder_list"
    failed_uploads="$log_folder/upload/fastq_failed_uploads.txt"
    files_to_upload="$log_folder/upload/fastq_files_to_upload.txt"
    touch "$failed_uploads" >/dev/null 2>&1
    touch "$files_to_upload" >/dev/null 2>&1
    > "$failed_uploads"

    for folder in $fastq_folder_list; do
        echo "Processing folder: $folder"
        sample=$(basename "$folder")
        > "$files_to_upload"
        sleep 5
        find "$folder" -mindepth 1 ! -type l \( -type f -o -type d \) > "$files_to_upload"

        while IFS= read -r file; do
            if [ ! -f "$file" ]; then
                echo "File does not exist: $file, skipping..."
                continue
            fi

            gs_path="gs://$bucket_name/fastq/$bcl/$sample/$(basename "$file")"
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
        done < "$files_to_upload"
    done

    if [ -s "$failed_uploads" ]; then
        echo -e "\nFailed to upload some files, please check $failed_uploads for details." >&2
    fi

    echo "Finished uploading fastq files."
else
    echo "No fastq folders found."
fi