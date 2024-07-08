# Snakefile
import os
import subprocess
import pandas as pd
from snakemake.utils import min_version

module_name = "slidetag_pipeline"
min_version("7.32.4")

# Read config parameters
def str_to_bool(s):
    return s.lower() in ["true", "1", "yes", "y"]

BASE_PATH = config['base_path']
SRC_PATH = config['src_path']
PYTHON_PATH = config['py_path']
BCL = config['bcl']
SHEET_ID = config['sheet_id']
REF_PATH = config['ref_path']
DEFAULT_BCL_PATH = config['bcl_path']
WORK_FOLDER = config['run_folder']

# Set run flags
RUN_MKFASTQ = str_to_bool(config['run_mkfastq'])
RUN_RNACOUNTS = str_to_bool(config['run_RNAcounts'])
RUN_CELLBENDER = str_to_bool(config['run_cellbender'])
RUN_SBCOUNTS = str_to_bool(config['run_SBcounts'])
RUN_SPATIAL = str_to_bool(config['run_spatial'])

# Set base path
BASE_FASTQ_PATH = os.path.join(BASE_PATH, BCL, 'fastq')
BASE_COUNT_PATH = os.path.join(BASE_PATH, BCL, 'count')
BASE_SPATIAL_PATH = os.path.join(BASE_PATH, BCL, 'spatial')
BASE_LOG_PATH = os.path.join(BASE_PATH, BCL, 'log')
CHECK_LOG = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'main/load_check.log')
SHEET_LOG = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'main/read_sheet.log')
BASE_MK_LOG_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'mkfastq_logs')
BASE_CR_LOG_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'counts_logs')
BASE_SB_LOG_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'spatial_logs')
paths_to_create = [BASE_MK_LOG_PATH, BASE_CR_LOG_PATH, BASE_SB_LOG_PATH]
for path in paths_to_create:
    if not os.path.exists(path):
        os.makedirs(path)


# run all
rule all:
    input:
        CHECK_LOG
    resources:
        mem_gb = config['read_check']['mem_gb'],
        runtime_min = config['read_check']['runtime_min']
    threads: config['read_check']['threads']

# read google sheets
rule read_sheet:
    log:
        log_file = SHEET_LOG
    params:
        bcl = BCL,
        sheet_id = SHEET_ID,
        data_path = BASE_PATH,
        python_path = PYTHON_PATH,
        src_path = SRC_PATH,
        bcl_path = DEFAULT_BCL_PATH,
        ref_path = REF_PATH,
        run_mkfastq = RUN_MKFASTQ,
        run_RNAcounts = RUN_RNACOUNTS,
        run_cellbender = RUN_CELLBENDER,
        run_SBcounts = RUN_SBCOUNTS,
        run_spatial = RUN_SPATIAL
    shell:
        """
        (
        set -euo pipefail
        umask 000
        src={params.src_path}/SubmitFuncs/read_info.py
        {params.python_path} $src {params.sheet_id} {params.bcl} {params.src_path} {params.data_path} {params.bcl_path} {params.ref_path} &> {log.log_file}
        sleep 5

        working_folder=$(echo "{log.log_file}" | sed 's|/main/read_sheet.log||g')
        if [[ {params.run_mkfastq} == 'True' ]]; then
            if [ -n "$(find "$working_folder/input/bcl2fastq/split_lanes" -type f >/dev/null 2>&1)" ]; then
                echo -e "\nERROR: $working_folder/input/bcl2fastq/split_lanes is empty." | tee -a {log.log_file}
            fi
            if [ -n "$(find "$working_folder/input/bcl2fastq/split_mkfastq" -type f >/dev/null 2>&1)" ]; then
                echo -e "\n$working_folder/input/split_mkfastq is empty." | tee -a {log.log_file}
            fi
        fi

        if [[ {params.run_RNAcounts} == 'True' || {params.run_cellbender} == 'True' ]]; then
            if [ -n "$(find "$working_folder/input/split_counts" -type f >/dev/null 2>&1)" ]; then
                echo -e "\nERROR: $working_folder/input/split_counts is empty." | tee -a {log.log_file}
            fi
        fi

        if [[ {params.run_SBcounts} == 'True' || {params.run_spatial} == 'True' ]]; then
            if [ -n "$(find "$working_folder/input/split_spatial" -type f >/dev/null 2>&1)" ]; then
                echo -e "\nERROR: $working_folder/input/split_spatial is empty." | tee -a {log.log_file}
            fi
        fi 
        ) || exit 1 || TRUE
        """

# check files
rule check_files:
    input:
        read_sheet = SHEET_LOG,
        base_fastq_path = BASE_FASTQ_PATH,
        base_count_path = BASE_COUNT_PATH,
        base_spatial_path = BASE_SPATIAL_PATH
    params:
        bcl = BCL,
        run_mkfastq = RUN_MKFASTQ,
        run_RNAcounts = RUN_RNACOUNTS,
        run_cellbender = RUN_CELLBENDER,
        run_SBcounts = RUN_SBCOUNTS,
        run_spatial = RUN_SPATIAL
    log:
        check_log = CHECK_LOG
    shell:
        """
        set -euo pipefail
        umask 000

        sheet_folder=$(awk '/Input sheets saved in:/ {{print $NF}}' {input.read_sheet})
        index_to_name="$sheet_folder/name_to_index.csv"
        bcl_path=$(awk '/BCL Path:/ {{print $NF}}' {input.read_sheet})
        sample_name_list=$(awk -F': ' '/RNA Indexes:/ {{print $2}}' {input.read_sheet} | sed 's/, / /g')

        if [[ {params.run_mkfastq} == 'True' ]] && [ ! -d "$bcl_path/{params.bcl}/Data/" ]; then
            echo "ERROR: running mkfastq but BCL data doesn't exist" | tee -a {log.check_log}
            exit 0
        fi

        for index in $sample_name_list; do
            echo "Checking $index" >> {log.check_log}
            fastq_path="{input.base_fastq_path}/$index"
            count_path="{input.base_count_path}/$index/outs/raw_feature_bc_matrix.h5"
            spatial_path="{input.base_spatial_path}/$index/SBcounts/SBcounts.h5"

            all_rna_fastq_exist=true
            all_sb_fastq_exist=true
            tail -n +2 "$index_to_name" | while IFS=',' read -r name rna_index vdj_index sb_index puck_id store_path; do
                if [[ ! "$rna_index" = "X" && ! "$rna_index" = "x" ]] && [[ -n "$rna_index" ]]; then
                    echo "Checking fastq files for $rna_index" &>> {log.check_log}
                    if [ ! -n "$(find "{input.base_fastq_path}/$rna_index" -type f >/dev/null 2>&1)" ]; then
                        all_rna_fastq_exist=false
                    fi
                fi
                if [[ ! "$sb_index" = "X" && ! "$sb_index" = "x" ]] && [[ -n "$sb_index" ]]; then
                    if [ ! -n "$(find ""{input.base_fastq_path}/$sb_index"" -type f >/dev/null 2>&1)" ]; then
                        all_sb_fastq_exist=false
                    fi
                fi
            done

            if [ ! -d "$count_path" ]; then
                counts_exist=false
            else
                counts_exist=true
            fi

            # RNAcounts checks
            if [[ {params.run_RNAcounts} == 'True' ]] && ([[ {params.run_mkfastq} == 'False' ]] && ! $all_rna_fastq_exist); then
                echo -e "\nERROR: running RNAcounts but fastq data doesn't exist for $index" | tee -a {log.check_log}
                continue
            fi

            # Cellbender checks
            if [[ {params.run_cellbender} == 'True' ]] && ([[ {params.run_RNAcounts} == 'False' ]] && ! $counts_exist); then
                echo -e "\nERROR: running cellbender but RNAcounts data doesn't exist for $index" | tee -a {log.check_log}
                continue
            fi

            # SBcounts checks
            if [[ {params.run_SBcounts} == 'True' ]] && ([[ {params.run_mkfastq} == 'False' ]] && ! $all_sb_fastq_exist); then
                echo -e "\nERROR: running SBcounts but SB fastqs doesn't exist for $index" | tee -a {log.check_log}
                continue
            fi

            # Spatial data checks
            if [[ {params.run_spatial} == 'True' ]] && ([[ {params.run_SBcounts} == 'False' ]] && [[ ! -f "$spatial_path" ]]); then
                echo -e "\nERROR: running spatial analysis but SBcounts.h5 doesn't exist for $index" | tee -a {log.check_log}
            fi
            
        done || TRUE
        """
