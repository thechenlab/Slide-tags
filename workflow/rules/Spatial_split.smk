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
CONDA_PATH = config['conda_path']
ENV_PATH = config['env_path']
REF_PATH = config['ref_path']
SRC_PATH = config['src_path']
BCL = config['bcl']
MATCH_RNA_LIST = config['rna_list']
WORK_FOLDER = config['run_folder']

# get file size
def get_file_size_in_mb(file_path):
    file_size_bytes = os.path.getsize(file_path)
    file_size_mb = file_size_bytes / (1024 * 1024)
    return file_size_mb

# Get resources based on file size
def get_resources(sample):
    path1 = os.path.join(BASE_COUNT_PATH, sample, 'outs/filtered_feature_bc_matrix.h5')
    path2 = os.path.join(BASE_COUNT_PATH, sample, 'cellbender_outs/cellbender_output_filtered.h5')
    if os.path.exists(path1):
        sample_path = path1
    else:
        sample_path = path2

    folder_size = get_file_size_in_mb(sample_path)
    if folder_size < 20:
        time = 3*60*60
    elif 30 < folder_size < 20:
        time = 4*60*60
    elif 50 > folder_size > 30:
        time = 5*60*60
    elif 70 > folder_size > 50:
        time = 6*60*60
    elif 90 > folder_size > 70:
        time = 7*60*60
    else:
        time = 10*60*60

    return time


# Set paths of log files
BASE_LOG_PATH = os.path.join(BASE_PATH, BCL, 'log')
BASE_COUNT_PATH = os.path.join(BASE_PATH, BCL, 'count')
BASE_SPATIAL_PATH = os.path.join(BASE_PATH, BCL, 'spatial')
SHEET_LOG = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'read_sheet.log')
BASE_SB_LOG_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'spatial_logs')
BASE_SHEET_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'input')
INDEX_TO_NAME = os.path.join(BASE_SHEET_PATH, 'name_to_index.csv') 

# Get list of log file
LOG_FILES = expand(os.path.join(BASE_SB_LOG_PATH, 'Spatial_{sample}.log'), sample=MATCH_RNA_LIST)
[os.remove(log) for log in LOG_FILES if os.path.exists(log)]

# Snakemake rule
rule all:
    input:
        LOG_FILES

rule spatial:
    input:
        name2_csv = INDEX_TO_NAME,
        spatial_path = BASE_SPATIAL_PATH
    log:
        spatial_log = os.path.join(BASE_SB_LOG_PATH, 'Spatial_{n}.log')
    params:
        bcl = BCL,
        sub_sample = lambda wildcards: wildcards.n,
        conda_path = CONDA_PATH,
        base_path = BASE_PATH,
        env_path = ENV_PATH,
        src_path = SRC_PATH,
        ref_path = REF_PATH
    resources:
        mem_gb = config['Spatial']['mem_gb'],
        disk_mb = config['Spatial']['disk_mb'],
        runtime_min = config['Spatial']['runtime_min']
    threads: config['Spatial']['threads']
    wildcard_constraints:
        n = '|'.join(MATCH_RNA_LIST)
    shell:
        """
        (
        set -euo pipefail
        ulimit -v unlimited -d unlimited
        umask 000
        
        export PATH="${params.conda_path}:$PATH"
        export PATH="{params.env_path}:$PATH"
        export PATH="{params.env_path}/bin:$PATH"
        export PATH="{params.env_path}/bin/R:$PATH"
        export R_LIBS="{params.env_path}/lib/R/library"
        export R_FUNC="{params.src_path}/Positioning"
        export REF_PATH="{params.ref_path}"
        export DATA_PATH="{params.base_path}"
        spatial_src="{params.src_path}/Positioning/run_spatial.R"

        while IFS=',' read -r name rna_index vdj_index sb_index puck_id store_path; do
            if [[ ! "$sb_index" = "X" && ! "$sb_index" = "x" ]] && [[ -n "$sb_index" ]]; then
                if [[ "$rna_index" == "{params.sub_sample}" ]]; then
                    target_file="{input.spatial_path}/$rna_index/Positions/seurat.qs"
                    if [ -f "$target_file" ]; then 
                        echo "seurat.qs for $rna_index already exists" | tee -a {log.spatial_log}
                        continue
                    else 
                        echo "Running spatial positioning for $rna_index" | tee -a {log.spatial_log}  
                        Rscript $spatial_src {params.bcl} "[$rna_index]" {threads} &>> {log.spatial_log}

                        sleep 10
                        if [ -f "$target_file" ]; then 
                            rm -rf "{input.spatial_path}/$rna_index/Positions/matrix.csv.gz" > /dev/null 2>&1
                            rm -rf "{input.spatial_path}/$rna_index/Positions/spatial_metadata.json" > /dev/null 2>&1
                            echo -e "SUCCESS" | tee -a {log.spatial_log}
                        else
                            echo -e "FAILURE" | tee -a {log.spatial_log}
                        fi
                    fi
                fi
            fi
        done < <(tail -n +2 {input.name2_csv}) 
        ) || TRUE
        """

