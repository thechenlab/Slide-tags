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
REF_PATH = config['ref_path']
SRC_PATH = config['src_path']
BCL = config['bcl']
WORK_FOLDER = config['run_folder']
SAMPLE_LIST = config['sample']

# get folder size
def get_folder_size(path):
    total_size = 0
    has_symlink = False
    for dirpath, dirnames, filenames in os.walk(path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            if os.path.islink(fp):
                has_symlink = True
                continue
            total_size += os.path.getsize(fp)
    if has_symlink:
        total_size *= 2
    return total_size / (1024**3)

# Get resources based on fastq size
def get_resources(sample):
    sample_name = re.sub(r'count_|vdj_', "", sample)
    sample_path = os.path.join(BASE_FASTQ_PATH, sample_name)
    folder_size = get_folder_size(sample_path)
    if folder_size < 10:
        time_count = 5*60*60
        time_vdj = 2*60*60
    elif 10 < folder_size < 20:
        time_count = 6*60*60
        time_vdj = 3*60*60
    elif 30 > folder_size > 40:
        time_count = 7*60*60
        time_vdj = 4*60*60
    elif 50 > folder_size > 40:
        time_count = 8*60*60
        time_vdj = 5*60*60
    elif 70 > folder_size > 60:
        time_count = 9*60*60
        time_vdj = 6*60*60
    else:
        time_count = 16*60*60
        time_vdj = 8*60*60
    if sample.startswith('count_'):
        return 16384, 8192, time_count, 68
    elif sample.startswith('vdj_'):
        return 16384, 8192, time_vdj, 68

# Set base path
BASE_FASTQ_PATH = os.path.join(BASE_PATH, BCL, 'fastq')
BASE_COUNT_PATH = os.path.join(BASE_PATH, BCL, 'count')
BASE_LOG_PATH = os.path.join(BASE_PATH, BCL, 'log')
SHEET_LOG = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'main/read_sheet.log')
BASE_CR_LOG_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'counts_logs')
BASE_SHEET_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'input')

# Get index file and log file
SAMPLE_INDEX = [n.replace("count_", "").replace("vdj_", "") for n in SAMPLE_LIST]
INDEX_FILES = expand(os.path.join(BASE_SHEET_PATH, "split_counts", "{n}.csv"), n=SAMPLE_LIST)
LOG_FILES = expand(os.path.join(BASE_CR_LOG_PATH, "{n}.log"), n=SAMPLE_LIST)
[os.remove(log) for log in LOG_FILES if os.path.exists(log)]


# Snakemake rule for RNAcounts 
rule all:
    input:
        LOG_FILES

rule RNAcounts:
    input:
        count_tsv = os.path.join(BASE_SHEET_PATH, "split_counts", "{n}.csv"),
        fastq_output = BASE_FASTQ_PATH,
        log_path = BASE_CR_LOG_PATH,
        count_path = BASE_COUNT_PATH,
        custom_ref_path = BASE_SHEET_PATH,
        read_sheet = SHEET_LOG
    log:
        RNAcounts_log = os.path.join(BASE_CR_LOG_PATH, "{n}.log")
    params:
        bcl = BCL,
        pkg_path = PKG_PATH,
        ref_path = REF_PATH,
        src_path = SRC_PATH
    resources:
        mem_gb = config['RNAcounts']['mem_gb'],
        disk_mb = config['RNAcounts']['disk_mb'],
        runtime_min = config['RNAcounts']['runtime_min']
    threads: config['RNAcounts']['threads']
    wildcard_constraints:
        n = '|'.join(SAMPLE_LIST)
    shell:
        """
        (
        set -euo pipefail
        ulimit -v unlimited -d unlimited
        umask 000

        export PATH="{params.pkg_path}/cellranger-arc-2.0.2/bin:$PATH"
        export PATH="{params.pkg_path}/cellranger-atac-2.1.0/bin:$PATH"
        
        cellranger_src={params.src_path}/Cellranger/run_cellranger.sh
        mkref_src={params.src_path}/Cellranger/make_refdata.sh
        config_path=$(echo {params.src_path} | sed 's/scripts/config/')/smk_config.yaml

        base_ref_path={params.ref_path}
        arc_ref="refdata-arc-GRCh38-2020-A, refdata-arc-mm10-2020-A"
        generate_bam=$(grep 'Generate BAM:' {input.read_sheet} | awk '{{print $NF}}')

        mkdir -p {input.count_path} 
        cd {input.count_path} 

        if [ ! -f {log.RNAcounts_log} ]; then
            touch {log.RNAcounts_log}
        fi

        while IFS=',' read -r idx ref_name chem vdj version; do
            idx=$(echo "$idx" | tr -d '\r')
            ref_name=$(echo "$ref_name" | tr -d '\r')
            dir_path="{input.count_path}/$idx"
            fastq_dir={input.fastq_output}/$idx
            lockfile="$dir_path/_lock"
            cmd_clean_files="rm -rf "$dir_path/SC_RNA_COUNTER_CS" \ 
                    "$dir_path/SC_ATAC_GEX_COUNTER_CS" \
                    "$dir_path/$idx.mri.tgz*" \
                    "$dir_path/_filelist*" \
                    "$dir_path/_finalstate*" \
                    "$dir_path/_mrosource*" \
                    "$dir_path/_jobmode*" \
                    "$dir_path/_perf*" \
                    "$dir_path/_sitecheck*" \
                    "$dir_path/_tags*" \
                    "$dir_path/_timestamp*" \
                    "$dir_path/_uuid*" \
                    "$dir_path/_vdrkill*" \
                    "$dir_path/_versions*" > /dev/null 2>&1"

            if [[ ! -d "$fastq_dir" ]]; then
                echo "Fastq folder $fastq_dir not found." | tee -a {log.RNAcounts_log}
                exit 0
            else
                echo "Fastq folder: $fastq_dir" &>> {log.RNAcounts_log}
                R1_files=$(ls $fastq_dir/$idx*_L00*_R1_*.fastq.gz 2> /dev/null || true)
                R2_files=$(ls $fastq_dir/$idx*_L00*_R2_*.fastq.gz 2> /dev/null || true)
                if [[ -z $R1_files || -z $R2_files ]]; then 
                    echo "Fastq files not found for $idx." | tee -a {log.RNAcounts_log}
                    exit 0
                fi
            fi

            if [[ -f "$dir_path/outs/web_summary.html" ]]; then
                echo "$idx analysis already completed, skipping..." | tee -a {log.RNAcounts_log}
                eval $cmd_clean_files
                exit 0
            fi

            if [[ ! -d "$dir_path/outs" ]] && [[ ! -d "$dir_path/tmp" ]] && [[ ! -f "{input.count_path}/__$idx.mro" ]]; then
                rm -rf "$dir_path" > /dev/null 2>&1
                sleep 30
                if [[ -d "$dir_path" ]]; then
                    echo "Failed to remove incomplete cellranger working folder due to permission issues." | tee -a {log.RNAcounts_log}
                    exit 0
                fi
            fi

            if [[ ! -d "$dir_path/outs" && ! -f "$lockfile" && ! "$dir_path/tmp" ]]; then
                rm -rf "$dir_path" > /dev/null 2>&1 || echo "Failed to remove incomplete cellranger working folder due to permission issues."
                rm -rf "*$idx*mro" > /dev/null 2>&1 || echo "Failed to remove incomplete cellranger working folder due to permission issues."
            fi

            if [[ -f "$lockfile" && -d "$dir_path/tmp" ]]; then
                rm -rf "$lockfile" > /dev/null 2>&1 || echo "Failed to delete cellranger lockfile due to permission issues."
            fi

            if echo "$arc_ref" | grep -q "$ref_name"; then
                atac_count="Y"
            else
                atac_count="N"
            fi

            echo "Index: $idx, Ref name: $ref_name" &>> {log.RNAcounts_log}
            reference="$base_ref_path/$ref_name"
            vdj_index="$vdj"

            custom_ref="{input.custom_ref_path}/custom_refdata.csv"
            if [[ -f $custom_ref ]]; then
                while IFS=',' read -r name index species gene; do
                    index=$(echo "$index" | tr -d '\r')
                    species=$(echo "$species" | tr -d '\r')
                    species=$(echo "$species" | tr '[:upper:]' '[:lower:]')
                    gene=$(echo "$gene" | tr -d '\r')
                    if [[ "$index" == "$idx" ]]; then
                        bash $mkref_src $gene $species &>> {log.RNAcounts_log}
                        created_ref="$base_ref_path/custom_refdata/${{gene}}/genome_${{species}}_${{gene}}"
                        if [[ -d "$created_ref" ]]; then
                            echo "Using custom reference with $gene." | tee -a {log.RNAcounts_log}
                            echo "Reference: $created_ref" &>> {log.RNAcounts_log}
                            reference="$created_ref"
                        else
                            echo "Reference data for $gene not found." | tee -a {log.RNAcounts_log}
                            exit 0
                        fi

                        break
                    fi
                done < $custom_ref
            fi

            if [[ "$version" == "V8" ]]; then
                export PATH="{params.pkg_path}/cellranger-8.0.1/bin:$PATH"
            else
                export PATH="{params.pkg_path}/cellranger-7.2.0/bin:$PATH"
            fi

            if [[ -n $chem ]]; then
                if [[ "$chem" == "ARC-v1" ]]; then
                    chemistry="ARC-v1" 
                else
                    chemistry="SC$chem"
                fi
            else
                chemistry="auto"
            fi

            log_file={log.RNAcounts_log}
            multiome_library="$(dirname {input.log_path})/input/split_counts/multiome_$idx.csv"
            count_path={input.count_path}
            bash $cellranger_src "$idx" "$reference" "$vdj_index" "$chemistry" "$version" "$log_file" "$fastq_dir" "$multiome_library" "$atac_count" "$dir_path" "$generate_bam" "$config_path"

            sleep 30
            if [[ -f "$dir_path/outs/web_summary.html" ]]; then
                eval $cmd_clean_files
                echo -e "SUCCESS" | tee -a $log_file
            else
                echo -e "FAILURE" | tee -a $log_file
            fi

        done < {input.count_tsv} 
        ) || TRUE
        """


