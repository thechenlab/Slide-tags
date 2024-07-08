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
PUCK_PATH = config['puck_path']
PUCK_IN_PATH = config['puck_in']
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
        total_size = 60*(1024**3)
    return total_size / (1024**3)

# Get resources based on fastq size
def get_resources(sample):
    sample_name = re.sub(r'SBcounts_', "", sample)
    sample_path = os.path.join(BASE_FASTQ_PATH, sample_name)
    folder_size = get_folder_size(sample_path)
    if folder_size < 10:
        time = 4*60*60
    elif 10 < folder_size < 20:
        time = 4*60*60
    elif 30 > folder_size > 40:
        time = 5*60*60
    elif 50 > folder_size > 40:
        time = 6*60*60
    elif 70 > folder_size > 60:
        time = 7*60*60
    else:
        time = 10*60*60
    return time


# Set base path
BASE_FASTQ_PATH = os.path.join(BASE_PATH, BCL, 'fastq')
BASE_SPATIAL_PATH = os.path.join(BASE_PATH, BCL, 'spatial')
BASE_LOG_PATH = os.path.join(BASE_PATH, BCL, 'log')
SHEET_LOG = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'read_sheet.log')
BASE_SB_LOG_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'spatial_logs')
BASE_SHEET_PATH = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'input')
INDEX_TO_NAME = os.path.join(BASE_SHEET_PATH, 'name_to_index.csv') 
SPATIAL = os.path.join(BASE_SHEET_PATH, 'Spatial.tsv')
SBCOUNTSV = os.path.join(BASE_SHEET_PATH, 'SBcounts.tsv') 

# Get list of log file and sample index
SAMPLE_LIST = [re.sub(r'SBcount_', "", n) for n in SAMPLE_LIST]
LOG_FILES = expand(os.path.join(BASE_SB_LOG_PATH, 'SBcount_{n}.log'), n=SAMPLE_LIST)
[os.remove(log) for log in LOG_FILES if os.path.exists(log)]

# Snakemake rule
rule all:
    input:
        LOG_FILES

rule SBcounts:
    input:
        puck_path = PUCK_PATH,
        puck_in = PUCK_IN_PATH,
        fastq_output = BASE_FASTQ_PATH,
        spatial_path = BASE_SPATIAL_PATH,
        spatial_tsv = SPATIAL,
        name2_csv = INDEX_TO_NAME
    log:
        SBcounts_log = os.path.join(BASE_SB_LOG_PATH, 'SBcount_{n}.log')
    params:
        sub_sample = lambda wildcards: wildcards.n,
        conda_path = CONDA_PATH,
        env_path = ENV_PATH,
        src_path = SRC_PATH,
        ref_path = REF_PATH
    resources:
        mem_gb = config['SBcounts']['mem_gb'],
        disk_mb = config['SBcounts']['disk_mb'],
        runtime_min = config['SBcounts']['runtime_min']
    threads: config['SBcounts']['threads']
    wildcard_constraints:
        n = '|'.join(SAMPLE_LIST)
    shell:
        """
        (
        set -euo pipefail
        ulimit -v unlimited -d unlimited
        umask 000

        export PATH="${params.conda_path}:$PATH"
        export PATH="{params.env_path}/bin:$PATH"
        export JULIA_PACKAGES_PATH="{params.env_path}/share/julia/packages"
        export JULIA_DEPOT_PATH="{params.env_path}/share/julia/tmp"
        export JULIA_PROJECT_PATH="{params.env_path}/julia/environments/env"
        export LIB_PUCK_PATH="{input.puck_path}"
        export LIB_PUCK_IN="{input.puck_in}"
        export JULIA_NUM_THREADS={threads}

        generate_puck_jl="{params.src_path}/SpatialCount/generate_puck_csv.jl"
        spatial_count_jl="{params.src_path}/SpatialCount/spatial_count.jl"

        cd {input.spatial_path} 
        touch {log.SBcounts_log}

        tail -n +2 {input.name2_csv} | while IFS=',' read -r name rna_index vdj_index sb_index puck_id store_path; do
            if [[ ! "$sb_index" = "X" && ! "$sb_index" = "x" ]] && [[ -n "$sb_index" ]]; then
                idx=$sb_index
                if [ $idx == {params.sub_sample} ]; then
                    working_dir="{input.spatial_path}/$rna_index/SBcounts"
                    mkdir -p $working_dir
                    cd $working_dir

                    if [[ -f $working_dir/SBcounts.h5 ]]; then 
                        echo "SBcounts.h5 for $rna_index exists, skipping..." &>> {log.SBcounts_log}
                        continue
                    fi

                    fastq_path="{input.fastq_output}/$idx"
                    puck_csv="{input.puck_path}/${{puck_id}}.csv" 
                    puck_path="$working_dir"

                    files_R3=("$fastq_path"/*_R3_*.fastq.gz)
                    if [ ${{#files_R3[@]}} -gt 0 ] && [ -e "${{files_R3[0]}}" ]; then
                        echo "Found _R3_.fastq.gz files, renaming..." &>> {log.SBcounts_log}
                        for f in "${{files_R3[@]}}"; do
                            f_r2=$(echo "$f" | sed 's/_R3_/_R2_/')
                            f_i2=$(echo "$f_r2" | sed 's/_R2_/_I2_/')
                            mv "$f_r2" "$f_i2" && echo "Renamed $f_r2 to $f_i2" &>> {log.SBcounts_log}
                            sleep 10
                            mv "$f" "$f_r2" && echo "Renamed $f to $f_r2" &>> {log.SBcounts_log}
                            sleep 10
                        done
                    fi

                    echo "Running SBcounts for $idx with $puck_id ..." | tee -a {log.SBcounts_log}
                    if [ ! -f $puck_csv ]; then
                        echo "$puck_id file doesn't exists, running generate_puck_csv.jl ..." &>> {log.SBcounts_log}
                        julia $generate_puck_jl &>> {log.SBcounts_log} || echo "Failed to run generate_puck_csv.jl" &>> {log.SBcounts_log}
                    fi

                    if [ ! -f $puck_csv ]; then
                        echo "$$puck_id file doesn't exists" &>> {log.SBcounts_log}
                        exit 0
                    else
                        echo "Copying puck csv file to spatial folder ..." &>> {log.SBcounts_log}
                        cp "$puck_csv" "$working_dir" &>/dev/null || echo "Failed to copy $puck_csv to $working_dir due to Permission denied" &>> {log.SBcounts_log}
                    fi

                    echo "Running spatial_count.jl ..." &>> {log.SBcounts_log}
                    julia $spatial_count_jl ${{fastq_path}} ${{working_dir}} &>> {log.SBcounts_log}

                    if [ -f $working_dir/SBcounts.h5 ]; then 
                        echo -e "SUCCESS" | tee -a {log.SBcounts_log}
                    else
                        echo -e "FAILURE" | tee -a {log.SBcounts_log}
                    fi
                fi
            fi
        done 
        ) || TRUE
        """
