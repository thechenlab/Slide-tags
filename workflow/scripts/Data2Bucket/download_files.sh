#!/bin/bash
umask 000

# Input parameters
bcl=$1
fastq_folder=$2
count_folder=$3
index_name=$4
pkg_path=$5
bucket_name=$6
index_fastq=$7

# Load packages
gcloud_pkg="$pkg_path/google-cloud-sdk/bin/gcloud"


### Download FASTQ files
if [ ! -z "$fastq_folder" ]; then
    echo -e "\n------------------------ Downloading FASTQ ------------------------ "
    gs_file_name="gs://$bucket_name/fastq/$bcl/$index_name"
    if "$gcloud_pkg" ls "$gs_file_name" 2>&1 | grep -q -v "CommandException: One or more URLs matched no objects."; then
        echo "$index_name fastqs exist in bucket, starting download..."
        mkdir -p "$fastq_folder/$index_name"

        download_fastqs() {
            "$gcloud_pkg" alpha storage cp -n "${gs_file_name}/$index_name$1" "$fastq_folder/$index_name/"
            wait $!
        }

        fastq_types=("_R1_*.fastq.gz" "_R2_*.fastq.gz")
        # Download I1 and I2 fastqs for cellranger count
        if [ "$index_fastq" == "true" ]; then
            fastq_types+=("_I1_*.fastq.gz" "_I2_*.fastq.gz")
            download_fastqs "*_I1_*.fastq.gz"
            download_fastqs "*_I2_*.fastq.gz"
        fi
        # Download R1 and R2 fastqs for spatial count
        download_fastqs "*_R1_*.fastq.gz"
        download_fastqs "*_R2_*.fastq.gz"

        # Check file sizes and re-download if necessary
        for fastq_type in "${fastq_types[@]}"; do
            local_files=("$fastq_folder/$index_name/"$fastq_type)
            for file in "${local_files[@]}"; do
                if [ -f "$file" ]; then
                    remote_size=$("$gcloud_pkg" du "${gs_file_name}/$file" | cut -f1)
                    local_size=$(stat -c%s "$file")
                    if [ "$remote_size" != "$local_size" ]; then
                        echo "Size mismatch for $file, re-downloading..."
                        download_fastqs "$fastq_type"
                    fi
                fi
            done
        done

    else
        echo "$index_name fastqs do not exist in bucket, skipping download..."
        exit 0
    fi
fi


### Download BAM files
if [ ! -z "$count_folder" ]; then
    echo -e "\n------------------------ Downloading BAM ------------------------ "
    gs_file_name="gs://$bucket_name/bam/$bcl/$index_name"
    if "$gcloud_pkg" ls "$gs_file_name" 2>&1 | grep -q -v "CommandException: One or more URLs matched no objects."; then
        echo "$index_name bams exist in bucket, starting download..."
        mkdir -p "$count_folder/$index_name/outs/"
        "$gcloud_pkg" alpha storage cp "$gs_file_name" "$count_folder/$index_name/outs/" > /dev/null 2>&1
    else
        echo "$index_name bams do not exist in bucket, skipping download..."
        exit 0
    fi
fi
