# Snakefile
import os
import re
import subprocess
import pandas as pd
from snakemake.utils import min_version

module_name = "slidetag_pipeline"
min_version("7.32.4")

# Read config parameters
BASE_PATH = config['base_path']
PKG_PATH = config['pkg_path']
BCL = config['bcl']
WORK_FOLDER = config['run_folder']
SAMPLE_LIST = config['sample']

# Get BCL store path for running mkfastq
def get_bcl_store_path(log_file_path):
    with open(log_file_path, 'r') as file:
        bcl_path = next((line.split()[2] for line in file if 'BCL Path' in line))
        return bcl_path

# Set base path
BASE_FASTQ_PATH = os.path.join(BASE_PATH, BCL, 'fastq')
BASE_LOG_PATH = os.path.join(BASE_PATH, BCL, 'log')
SHEET_LOG = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'main/read_sheet.log')
BASE_MK_LOG_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'mkfastq_logs')
BASE_SHEET_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'input')

# Get index file and log file
LOG_FILES = expand(os.path.join(BASE_MK_LOG_PATH, "{n}.log"), n=SAMPLE_LIST)
INDEX_LIST = expand(os.path.join(BASE_SHEET_PATH, "split_mkfastq/{n}.csv"), n=SAMPLE_LIST)
[os.remove(log) for log in LOG_FILES if os.path.exists(log)]


# Snakemake rule for RNA and ATAC mkfastq
rule all:
    input:
        LOG_FILES

rule mkfastq:
    input:
        Indexes =  os.path.join(BASE_SHEET_PATH, "split_mkfastq/{n}.csv"), 
        bcl2fq_Indexes = os.path.join(BASE_SHEET_PATH, "bcl2fastq/split_lanes/{n}.csv"), 
        fastq_output = BASE_FASTQ_PATH
    log:
        mkfastq_log = os.path.join(BASE_MK_LOG_PATH, "{n}.log")
    params:
        bcl = BCL,
        pkg_path = PKG_PATH,
        bcl_path = lambda wildcards: get_bcl_store_path(SHEET_LOG)
    resources:
        mem_gb = config['mkfastqt']['mem_gb'],
        disk_mb = config['mkfastqt']['disk_mb'],
        runtime_min = config['mkfastqt']['runtime_min']
    threads: config['mkfastqt']['threads']
    wildcard_constraints:
        n = '|'.join(SAMPLE_LIST)
    shell:
        """
        set -euo pipefail
        umask 000
        ulimit -v unlimited -d unlimited
        export PATH="{params.pkg_path}/bcl2fastq2_v2.20.0/bin:$PATH"
        
        job_type=$(basename {input.Indexes} | cut -d '.' -f 1)
        index_list=$(awk -F ',' 'NR>1 && $2 != "" {{print $2}}' {input.Indexes} | tr '\n' ' ')
        lane_num=$(echo "{input.Indexes}" | grep -o 'lane[0-9]*' | cut -d 'e' -f 2)

        if [[ "$job_type" == "ATAC_lane"* ]]; then
            lane="${{lane_num}}_atac" 
        elif [[ "$job_type" == "ATAC_Index_lane"* ]]; then
            lane="${{lanlane_nume}}_arc"
        else
            lane="${{lane_num}}"
        fi

        all_files_exist=true
        for index in $index_list; do
            R1_files=$(ls {input.fastq_output}/$index/$index*_L00${{lane_num}}_R1_*.fastq.gz 2> /dev/null || true)
            R2_files=$(ls {input.fastq_output}/$index/$index*_L00${{lane_num}}_R2_*.fastq.gz 2> /dev/null || true)
            if [[ -z $R1_files || -z $R2_files ]]; then 
                all_files_exist=false
                break
            fi
        done

        if [ "$all_files_exist" = false ]; then
            echo "Checking fastq files in mkfastq folder." &>> {log.mkfastq_log}
            target_path="{input.fastq_output}/mkfastq_$lane/outs/fastq_path"
            bcl2fq_stats="{input.fastq_output}/mkfastq_$lane/Stats/Stats.json"
            echo "Checking $target_path" &>> {log.mkfastq_log}

            fastq_folder="" 
            if [ -d "$target_path" ]; then
                fastq_folder=$(find "$target_path" -maxdepth 1 -type d | grep -vE "(Reports|Stats)" | head -3 | tail -1)
            elif [ -f "$bcl2fq_stats" ]; then
                fastq_folder="{input.fastq_output}/mkfastq_$lane"
            fi

            if [ -z "$fastq_folder" ]; then
                for index in $index_list; do
                    if [ -d "$fastq_folder" ] && [ -f "$fastq_folder/$index*.fastq.gz" ]; then
                        echo "Found fastq files in mkfastq folder for $index. Proceeding with file moving." | tee -a {log.mkfastq_log}
                        for file in $(find "$fastq_folder" -type f -name "*.fastq.gz"); do
                            mv "$file" "{input.fastq_output}"
                        done
                        mkdir -p {input.fastq_output}/$index
                        mv {input.fastq_output}/$index*.fastq.gz {input.fastq_output}/$index
                    fi
                done
            fi
        else
            echo "All required fastq.gz files exist. Skipping mkfastq processing." | tee -a {log.mkfastq_log}
            exit 0
        fi

        if [ -d {input.fastq_output}/mkfastq_$lane/tmp ] ; then
            rm -rf {input.fastq_output}/mkfastq_$lane
        fi
        mkdir -p {input.fastq_output}
        cd {input.fastq_output}

        echo "Running cellranger mkfastq" | tee -a {log.mkfastq_log}
        if [[ "$job_type" == "ATAC_Index_lane"* ]]; then
            export PATH="{params.pkg_path}/cellranger-arc-2.0.2/bin:$PATH"
            mkfastq_cmd="cellranger-arc mkfastq"
            mask_info="Y50,I8,Y24,Y50"
        elif [[ "$job_type" == "ATAC_lane"* ]]; then
            export PATH="{params.pkg_path}/cellranger-arc-2.0.2/bin:$PATH"
            mkfastq_cmd="cellranger-arc mkfastq"
            mask_info="Y50,I8,Y24,Y50"
        else
            export PATH="{params.pkg_path}/cellranger-8.0.1/bin:$PATH"
            mkfastq_cmd="cellranger mkfastq"
            mask_info="Y28n*,I10,I10,Y90n*"
        fi

        (time stdbuf -oL -eL $mkfastq_cmd \
            --run={params.bcl_path}{params.bcl} \
            --id="mkfastq_$lane" --csv={input.Indexes} \
            --delete-undetermined &>> {log.mkfastq_log}) &>> {log.mkfastq_log} || true

        if [ ! -d "{input.fastq_output}/mkfastq_$lane/outs" ]; then
            rm -rf {input.fastq_output}/mkfastq_$lane &>> {log.mkfastq_log}
            echo "$mkfastq_cmd failed, switching to run bcl2fastq" | tee -a {log.mkfastq_log}

            bcl_run_xml="{params.bcl_path}{params.bcl}/RunInfo.xml"
            mask_num=$(xmllint --xpath '//Read' "$bcl_run_xml" 2>/dev/null | grep -oP 'Read Number="\d+" NumCycles="\d+"' | awk -F'"' '{{print $4}}' | paste -sd ',' -)

            if [[ "$mask_num" == "170,10,24,134" ]]; then
                mask_info="Y170,I8n*,Y24,Y134"
            fi

            if [[ "$mask_num" == "150,10,24,150" ]]; then
                if [[ "$job_type" == "ATAC_"* ]]; then
                    mask_info="Y150,I8n*,Y24,Y150"
                else
                    mask_info="Y150,I10,Y24,Y150"
                fi
            fi

            echo "Using bases-mask: $mask_info" &>> {log.mkfastq_log}
            if [[ -f {input.bcl2fq_Indexes} ]]; then
                echo "Using bcl2fastq Indexes file" &>> {log.mkfastq_log}
            else
                echo "bcl2fastq Indexes file not found" | tee -a {log.mkfastq_log}
                exit 1
            fi
            
             (time stdbuf -oL -eL bcl2fastq \
                --use-bases-mask=$mask_info \
                --create-fastq-for-index-reads \
                --minimum-trimmed-read-length=8 \
                --mask-short-adapter-reads=8 \
                --ignore-missing-positions \
                --ignore-missing-controls \
                --ignore-missing-filter \
                --ignore-missing-bcls \
                --tiles $lane_num \
                -R {params.bcl_path}{params.bcl} \
                --output-dir={input.fastq_output}/mkfastq_$lane \
                --sample-sheet={input.bcl2fq_Indexes} &>> {log.mkfastq_log} ) &>> {log.mkfastq_log} || true
            if [ ! -f "{input.fastq_output}/mkfastq_$lane/Stats/Stats.json" ]; then
                echo "bcl2fastq mkfastq also failed" | tee -a {log.mkfastq_log}

            fi
        fi

        sleep 30
        if [ -f {input.fastq_output}/mkfastq_$lane/outs/fastq_path/Reports/html/index.html ]; then
            echo "Moving fastq files from cellranger output" &>> {log.mkfastq_log}

            rm -rf "{input.fastq_output}/mkfastq_$lane/MAKE_FASTQS_CS" &>> {log.mkfastq_log}
            rm -rf "{input.fastq_output}/mkfastq_$lane/mkfastq_$lane.mri.tgz" &>> {log.mkfastq_log}
            rm -rf "{input.fastq_output}/mkfastq_$lane/outs/interop_path" &>> {log.mkfastq_log}
            
            target_path="{input.fastq_output}/mkfastq_$lane/outs/fastq_path"
            fastq_folder=$(find "$target_path" -maxdepth 1 -type d | grep -vE "(Reports|Stats)"  | head -3 | tail -1)
            if [ -d $fastq_folder ]; then
                max_attempts=3
                for index in $index_list; do
                    echo "Moving fastq files for $index" &>> {log.mkfastq_log}
                    mkdir -p "{input.fastq_output}/$index"
                    mv_attempts=0
                    if ls "$fastq_folder"/"$index"* 1> /dev/null 2>&1; then
                        mv_attempts=0
                        until mv "$fastq_folder"/"$index"* "{input.fastq_output}/$index" || [ $mv_attempts -ge $max_attempts ]; do
                            mv_attempts=$((mv_attempts + 1))
                            echo "Attempt $mv_attempts to move files for $index failed, retrying..." &>> {log.mkfastq_log}
                            sleep 5
                        done
                    fi
                    if [ $mv_attempts -ge $max_attempts ]; then
                        echo "Failed to move files for $index after $mv_attempts attempts" &>> {log.mkfastq_log}
                    else
                        echo "SUCCESS" | tee -a {log.mkfastq_log}
                    fi
                done
            fi
        elif [ -f {input.fastq_output}/mkfastq_$lane/Stats/Stats.json ]; then
            echo "Moving fastq files from bcl2fastq output" &>> {log.mkfastq_log}
            fastq_folder="{input.fastq_output}/mkfastq_$lane"
            rm -rf "$fastq_folder"/"Undetermined"* &>> {log.mkfastq_log}

            for index in $index_list; do
                mkdir -p "{input.fastq_output}/$index"
                max_attempts=3
                if ls "$fastq_folder"/"$index"* 1> /dev/null 2>&1; then
                    mv_attempts=0
                    until mv "$fastq_folder"/"$index"* "{input.fastq_output}/$index" || [ $mv_attempts -ge $max_attempts ]; do
                        mv_attempts=$((mv_attempts + 1))
                        echo "Attempt $mv_attempts to move files for $index failed, retrying..." &>> {log.mkfastq_log}
                        sleep 5
                    done
                fi
                if [ $mv_attempts -ge $max_attempts ]; then
                    echo "Failed to move files for $index after $mv_attempts attempts" &>> {log.mkfastq_log}
                else
                    echo -e "SUCCESS" | tee -a {log.mkfastq_log}
                fi
            done

            rm -rf "{input.fastq_output}/Undetermined"
        else
            echo "Failed to generate fastq files" &>> {log.mkfastq_log}
            echo -e "FAILURE" | tee -a {log.mkfastq_log}
        fi || TRUE
        """