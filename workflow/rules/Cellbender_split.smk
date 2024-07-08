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
ENV_PATH = config['env_path']
BCL = config['bcl']
WORK_FOLDER = config['run_folder']
SAMPLE_LIST = config['sample']

# get sample size from cellranger outs
def get_file_size_in_mb(file_path):
    file_size_bytes = os.path.getsize(file_path)
    file_size_mb = file_size_bytes / (1024 * 1024)
    return file_size_mb

# Get resources based on raw h5 file size
def get_resources(sample):
    sample_name = re.sub(r'count_|multiome_', "", sample)
    sample_path = os.path.join(BASE_COUNT_PATH, sample_name, 'outs/raw_feature_bc_matrix.h5')
    folder_size = get_file_size_in_mb(sample_path)
    if folder_size < 10:
        time = 12*60*60
    elif 10 < folder_size < 30:
        time = 14*60*60
    elif 30 > folder_size > 50:
        time = 16*60*60
    elif 50 > folder_size > 70:
        time = 18*60*60
    elif 70 > folder_size > 90:
        time = 20*60*60
    elif 90 > folder_size > 110:
        time = 22*60*60
    elif 110 > folder_size > 130:
        time = 24*60*60
    elif 130 > folder_size > 150:
        time = 26*60*60
    else:
        time = 30*60*60
    return time

# set base path
BASE_COUNT_PATH = os.path.join(BASE_PATH, BCL, 'count')
BASE_LOG_PATH = os.path.join(BASE_PATH, BCL, 'log')
SHEET_LOG = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'main/read_sheet.log')
BASE_CR_LOG_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'counts_logs')
BASE_SHEET_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'input')

# Get list of log file
SAMPLE_LIST = [n.replace("count_", "").replace("multiome_", "").replace('.csv', '') for n in SAMPLE_LIST]
LOG_FILES = expand(os.path.join(BASE_CR_LOG_PATH, "Cellbender_{n}.log"), n=SAMPLE_LIST)
[os.remove(log) for log in LOG_FILES if os.path.exists(log)]


# Snakemake rule for each cellbender 
rule all:
    input:
        LOG_FILES

rule cellbender:
    input:
        count_path = BASE_COUNT_PATH,
    log:
        cellbender_log = os.path.join(BASE_CR_LOG_PATH, "Cellbender_{n}.log")
    params:
        env_path = ENV_PATH,
        expected_cells = config["expected_cells"],
        total_droplets = config["total_droplets_included"],
        sub_sample =  lambda wildcards: wildcards.n
    resources:
        mem_gb = config['Cellbender']['mem_gb'],
        disk_mb = config['Cellbender']['disk_mb'],
        runtime_min = config['Cellbender']['runtime_min']
    threads: config['Cellbender']['threads']
    wildcard_constraints:
        n = '|'.join(SAMPLE_LIST)
    shell:
        """
        set -euo pipefail
        umask 000
        ulimit -v unlimited -d unlimited
        export PATH="$(echo {params.env_path} | sed 's/main/cellbender2/')/bin/:$PATH"

        cd {input.count_path} 
        idx="{params.sub_sample}"
        echo "Processing $idx" &>> {log.cellbender_log}
        cellranger_file="{input.count_path}/$idx/outs/raw_feature_bc_matrix.h5"
        cellranger_bc="{input.count_path}/$idx/outs/raw_feature_bc_matrix/barcodes.tsv.gz"
        cellbender_path="{input.count_path}/$idx/cellbender_outs"

        if [ -f "$cellbender_path/cellbender_output.h5" ]; then
            echo "Cellbender output already exists for $idx" | tee -a {log.cellbender_log}
            exit 0
        else
            rm -rf $cellbender_path
        fi

        if [ -f "$cellranger_file" ]; then
            echo "Running cellbender for $idx" | tee -a {log.cellbender_log}
            cmd="time stdbuf -oL -eL cellbender remove-background --input $cellranger_file --output $cellbender_path/cellbender_output.h5"

            # add expected-cells
            if [ "{params.expected_cells}" != "none" ] && [ "{params.expected_cells}" != "None" ]; then 
                cmd+=" --expected-cells {params.expected_cells}"
            fi
            
            # add total-droplets-included
            raw_cell_num=$(wc -l < "$cellranger_bc" | xargs)
            if [ "{params.total_droplets}" != "none" ] && [ "{params.total_droplets}" != "None" ]; then
                cmd+=" --total-droplets-included {params.total_droplets}"
            elif [ "$raw_cell_num" -gt 30000 ]; then
                cmd+=" --total-droplets-included 50000"
            elif [ "$raw_cell_num" -gt 20000 ]; then
                cmd+=" --total-droplets-included 40000"
            fi
            
            echo "Running cellbender for $idx with command: $cmd" >> {log.cellbender_log}
            mkdir -p {input.count_path}/$idx/cellbender_outs 

            # run cellbender
            cd $cellbender_path
            echo "$cmd" >> {log.cellbender_log}
            eval $cmd &>> {log.cellbender_log} || echo "cellbender failed" &>> {log.cellbender_log}

            # check if cellbender output exists
            if [ -f "$cellbender_path/cellbender_output.h5" ]; then
                echo -e "SUCCESS" | tee -a {log.cellbender_log}
            else
                echo -e "FAILURE" | tee -a {log.cellbender_log}
            fi
        fi || TRUE
        """