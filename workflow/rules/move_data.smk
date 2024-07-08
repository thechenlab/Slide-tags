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
BCL = config['bcl']
WORK_FOLDER = config['run_folder']

# Set paths of log files
BASE_COUNT_PATH = os.path.join(BASE_PATH, BCL, 'count')
BASE_SPATIAL_PATH = os.path.join(BASE_PATH, BCL, 'spatial')
BASE_LOG_PATH = os.path.join(BASE_PATH, BCL, 'log')
MOVE_LOG = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'main/move_data.log')
INDEX_TO_NAME = os.path.join(BASE_LOG_PATH, WORK_FOLDER, 'input/name_to_index.csv')

# Snakemake rule
rule move_data:
    input:
        spatial_path = BASE_SPATIAL_PATH,
        count_path = BASE_COUNT_PATH,
        name2_csv = INDEX_TO_NAME,
    log:
        move_log = MOVE_LOG
    resources:
        mem_gb = config['move_data']['mem_gb'],
        disk_mb = config['move_data']['disk_mb'],
        runtime_min = config['move_data']['runtime_min']
    threads: config['move_data']['threads']
    shell:
        """
        set -euo pipefail
        umask 000

        touch {log.move_log}
        
        # Move spatial data
        while IFS=',' read -r name rna_index vdj_index sb_index puck_id store_path; do

            if [[ -z "$store_path" ]] && [[ ! "$sb_index" = "X" && ! "$sb_index" = "x" ]]; then
                echo "Store path for $name is empty. Skipping..." >> {log.move_log}
                continue
            fi
            
            if [[ "$store_path" == */ ]]; then
                final_folder="$store_path$name"
            else
                final_folder="$store_path/$name"
            fi

            echo -e "\nCreated directory: $final_folder" &>> {log.move_log}
            mkdir -p "$final_folder" 2>&1 || echo -e "Failed to create directory $final_folder. \nCheck permissions and path correctness." &>> {log.move_log}
            
            # files to move 
            counts_file="{input.count_path}/$rna_index/outs/filtered_feature_bc_matrix.h5"
            web_summary="{input.count_path}/$rna_index/outs/web_summary.html"
            vdj_file="{input.count_path}/$rna_index/outs/clonotypes.csv"
            vdj_contig="{input.count_path}/$rna_index/outs/all_contig_annotations.csv"
            vdj_consensus="{input.count_path}/$rna_index/outs/consensus_annotations.csv"
            vdj_filter="{input.count_path}/$rna_index/outs/filtered_contig_annotations.csv"
            vdj_metric="{input.count_path}/$rna_index/outs/metrics_summary.csv"
            cellbender_file="{input.count_path}/$rna_index/cellbender_outs/cellbender_output_filtered.h5"
            spatial_folder="{input.spatial_path}/$rna_index/Positions"

            if [[ -f "$counts_file" ]]; then
                cp "$counts_file" "$final_folder" 2>&1 || echo -e "Failed to copy filtered_feature_bc_matrix.h5. \nCheck permissions and path correctness." &>> {log.move_log}
                echo "Copied $counts_file to created folder." &>> {log.move_log}
                cp "$web_summary" "$final_folder" 2>&1 || echo -e "Failed to copy counts_web_summary. \nCheck permissions and path correctness." &>> {log.move_log}
                echo "Copied $web_summary to created folder." &>> {log.move_log}
            else
                echo "Counts file not found for $name. Skipping..." &>> {log.move_log}
            fi

            if [[ -f "$vdj_file" ]]; then
                cp "$vdj_file" "$final_folder" 2>&1 || echo -e "Failed to copy clonotypes.csv. \nCheck permissions and path correctness." &>> {log.move_log}
                echo "Copied $vdj_file to created folder." &>> {log.move_log}
                cp "$web_summary" "$final_folder" 2>&1 || echo -e "Failed to copy vdj_web_summary. \nCheck permissions and path correctness." &>> {log.move_log}
                echo "Copied $web_summary to created folder." &>> {log.move_log}
                cp "$vdj_contig" "$final_folder" 2>&1 || echo -e "Failed to copy all_contig_annotations.csv. \nCheck permissions and path correctness." &>> {log.move_log}
                echo "Copied $vdj_contig to created folder." &>> {log.move_log}
                cp "$vdj_consensus" "$final_folder" 2>&1 || echo -e "Failed to copy consensus_annotations.csv. \nCheck permissions and path correctness." &>> {log.move_log}
                echo "Copied $vdj_consensus to created folder." &>> {log.move_log}
                cp "$vdj_filter" "$final_folder" 2>&1 || echo -e "Failed to copy filtered_contig_annotations.csv. vCheck permissions and path correctness." &>> {log.move_log}
                echo "Copied $vdj_filter to created folder." &>> {log.move_log}
                cp "$vdj_metric" "$final_folder" 2>&1 || echo -e "Failed to copy metrics_summary.csv. \nCheck permissions and path correctness." &>> {log.move_log}
                echo "Copied $vdj_metric to created folder." &>> {log.move_log}
            else
                echo "VDJ file not found for $name. Skipping..." &>> {log.move_log}
            fi

            if [[ -f "$cellbender_file" ]]; then
                cp "$cellbender_file" "$final_folder" 2>&1 || echo -e "Failed to copy cellbender_output_filtered.h5. \nCheck permissions and path correctness." &>> {log.move_log}
                echo "Copied $cellbender_file to created folder." &>> {log.move_log}
            else
                echo "Cellbender file not found for $name. Skipping..." &>> {log.move_log}
            fi

            if [[ -d "$spatial_folder" ]]; then
                cp -r "$spatial_folder" "$final_folder" 2>&1 || echo -e "Failed to copy spatial results. \nCheck permissions and path correctness." &>> {log.move_log}
                echo "Copied $spatial_folder outputs to created folder." &>> {log.move_log}
            else
                echo "Spatial output not found for $name. Skipping..." &>> {log.move_log}
            fi

        done < <(tail -n +2 {input.name2_csv}) || TRUE
        """
