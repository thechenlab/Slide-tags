#!/bin/bash
umask 000
current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source "$current_dir/config/config.sh"
source activate "$ENV_PATH"


# Function to print help information
print_help() {
    echo "Usage: Slidetag Pipeline [options]"
    echo ""
    echo "Required:"
    echo "  -bcl [value]                  Set BCL name."
    echo "Options:"
    echo "  -h                            Display this help message."
    echo "  -ra, -run_all                 Run mkfastq, RNAcounts, SBcounts and Spatial analysis."
    echo "  -mk, -run_mkfastq             Run cellranger mkfastq or bcl2fastq."
    echo "  -cr, -run_RNAcounts           Run cellranger count or cellranger arc."
    echo "  -cb, -run_cellbender          Run cellbender based on Cellranegr count results."
    echo "  -sb, -run_SBcounts            Run Spatial beads counts."
    echo "  -sp, -run_spatial             Run Spatial analysis for cell positionings."
    echo "  -us, -use_sheet               Get input sheets form the current working run."
    echo "  -mv, -mv_file                 Move results to store path."
    echo "  -gb, -generate_bam            Generate bams when running Cellranger."
    echo "  -uf, -upload_fastq            Upload fastqs to google bucket."
    echo "  -ub, -upload_bam              Upload bams to google bucket."
    echo "  -rf, -rm_fastq                Remove local fastqs."
    echo "  -rb, -rm_bam                  Remove local bams."
    echo "  -df, -download_fastq          Download fastqs from google bucket."
    echo "  -db, -download_bam            Download bams from google bucket."
    echo "  -f, -force                    Force to re-run selected jobs."
    echo "  -ec [value], -expected_cells  From cellbender parameters."
    echo "  -td [value], -total_droplets_included"
    echo "FALSE or NONE at default for the above parameters."
    echo ""
}


# Initialize variables
use_sheet=false
run_all=false
run_mkfastq=false
run_RNAcounts=false
run_cellbender=false
run_SBcounts=false
run_spatial=false
mv_file=false
generate_bam=false
upload_fastq=false
upload_bam=false
download_fastq=false
download_bam=false
rm_fastq=false
rm_bam=false
force=false
expected_cells=""
total_droplets_included=""

# Parse input parameters
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -bcl) bcl="$2"; shift 2 ;;
        -ra|-run_all) run_all=true; shift ;;
        -mk|-run_mkfastq) run_mkfastq=true; shift ;;
        -cr|-run_RNAcounts) run_RNAcounts=true; shift ;;
        -cb|-run_cellbender) run_cellbender=true; shift ;;
        -sb|-run_SBcounts) run_SBcounts=true; shift ;;
        -sp|-run_spatial) run_spatial=true; shift ;;
        -us|-use_sheet) use_sheet=true; shift ;;
        -mv|-mv_file) mv_file=true; shift ;;
        -gb|-generate_bam) generate_bam=true; shift ;;
        -uf|-upload_fastq) upload_fastq=true; shift ;;
        -ub|-upload_bam) upload_bam=true; shift ;;
        -df|-download_fastq) download_fastq=true; shift ;;
        -db|-download_bam) download_bam=true; shift ;;
        -rf|-rm_fastq) rm_fastq=true; shift ;;
        -rb|-rm_bam) rm_bam=true; shift ;;
        -ec|-expected_cells) expected_cells="$2"; shift 2 ;;
        -td|-total_droplets_included) total_droplets_included="$2"; shift 2 ;;
        -f|-force) force=true; shift ;;
        -h) print_help; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
done


# check if bcl param is empty
if [ -z "$bcl" ]; then
    echo "Error: param for -bcl is empty"
    exit 1
else 
    echo -e "\n=================================="
    echo "    Running Slide-tags Pipeline   "
    echo "=================================="
    echo -e "\nProcessing BCL: $bcl"
    mkdir -p "$BASE_DATA_PATH/$bcl"
fi

# Set base paths and create directories
python_path="$ENV_PATH/bin/python"
julia_path="$ENV_PATH/bin/julia"
smk_config="$current_dir/config/smk_config.yaml"
base_fastq_path="$BASE_DATA_PATH/$bcl/fastq"
base_count_path="$BASE_DATA_PATH/$bcl/count"
base_spatial_path="$BASE_DATA_PATH/$bcl/spatial"
base_log_path="$BASE_DATA_PATH/$bcl/log"
base_smk_path="$WORKFLOW_PATH/workflow/rules"
base_src_path="$WORKFLOW_PATH/workflow/scripts"
source "$base_src_path/SubmitFuncs/submit_functions.sh"
mkdir -p "$base_fastq_path"
mkdir -p "$base_count_path"
mkdir -p "$base_spatial_path"
mkdir -p "$base_log_path"

# Get chunk size for each job
size_mk="$CHUNK_SIZE_MKFASTQ"
size_cr="$CHUNK_SIZE_RNACOUNTS"
size_cb="$CHUNK_SIZE_CELLBENDER"
size_sb="$CHUNK_SIZE_SBCOUNT"
size_sp="$CHUNK_SIZE_POSITION"


# check if run_all is true
if [[ "$run_all" = "true" ]]; then
    run_mkfastq=true
    run_RNAcounts=true
    run_SBcounts=true
    run_spatial=true
fi
# Check if the cluster path is valid
if [[ -z "$CLUSTER_PATH" ]] || [[ ! -e "$CLUSTER_PATH" ]]; then
    use_cluster=false
else
    use_cluster=true
fi


###### Upload fastq and bams to gcloud
if [ "$upload_fastq" = "true" ]; then
    upload_fastq "$bcl" "$base_fastq_path" "$base_log_path" "$PKG_PATH" "$GOOGLE_CLOUD_BUCKET" "$use_cluster"
    exit 0
fi
if [ "$upload_bam" = "true" ]; then
    upload_bam "$bcl" "$base_count_path" "$base_log_path" "$PKG_PATH" "$GOOGLE_CLOUD_BUCKET" "$use_cluster"
    exit 0
fi


###### Delete fastq and bams
if [ "$rm_fastq" = "true" ]; then
    echo -e "\n------------------------ Removing local FASTQ ------------------------ "
    rm -rf "$base_fastq_path" >/dev/null 2>&1 || echo "Failed to remove fastq due to permission issues."
    echo "Done."
    exit 0
fi
if [ "$rm_bam" = "true" ]; then
    echo -e "\n------------------------ Removing local BAMS ------------------------ "
    find "$base_count_path" -type f \( -name 'possorted_genome_bam.bam' -o -name 'possorted_genome_bam.bam.bai' \) -delete 2>/dev/null || echo "Failed to remove bam due to permission issues."
    find "$base_count_path" -type d -name 'SC_RNA_COUNTER_CS' -exec rm -rf {} + 2>/dev/null || echo "Failed to remove SC_RNA_COUNTER_CS due to permission issues."
    echo "Done."
    exit 0
fi


###### Set working directory
working_file="${base_log_path}/working"
if [ "$use_sheet" = "true" ]; then
    if [ -f "$working_file" ]; then
        log_folder=$(cat "$working_file")
    fi
else
    for ((i=0; i<100; i++)); do
        folder_name="${i}run"
        folder_to_create="${base_log_path}/${folder_name}"
        if [ ! -d "$folder_to_create" ]; then
            mkdir -p "$folder_to_create"
            if [ -f "$working_file" ]; then
                rm "$working_file" >/dev/null 2>&1
            fi
            touch "$working_file" >/dev/null 2>&1
            echo "$folder_to_create" > "$working_file"
            log_folder="$folder_to_create"
            break
        fi
    done
fi
echo -e "Current working folder: $log_folder"
mkdir -p "$log_folder/main" 
base_sheet_folder="$log_folder/input"
source "$base_src_path/SubmitFuncs/submit_functions.sh"


###### Read google sheet and check files
if [ "$use_sheet" = "false" ]; then
    echo -e "\n------------------------ Reading google sheet and checking files ------------------------ "
    submit_read_sheet_job "$use_cluster"
    cleanup_snakemake_logs "$(pwd)" "$log_folder" "$base_log_path"
    chmod -R 777 "$log_folder" >/dev/null 2>&1 
fi 


#---------------------------------------------------------------------------------------------------------------------------------------
###### Submit Cellranger mkfastq jobs
if [ "$run_mkfastq" = "true" ]; then
    index_list=$(make_mkfastq_input_list "$base_fastq_path" "$log_folder")
    index_list_check=($(echo "${index_list[@]}" | grep -v '^\s*$'))
    if [ "${#index_list_check[@]}" -gt 0 ] || [ "${#index_list_check[@]}" -eq 0 ]; then
        echo -e "\n------------------------  Running cellranger mkfastq ------------------------ "
        check_list=()
        for item in "${index_list_check[@]}"; do
            check_list+=("$(echo "$item" | sed -e 's/^\(ATAC_\)\?Index_lane/mkfastq_/')")
        done
        submit_batch_jobs "$use_cluster" "$index_list" "$size_mk" "submit_mkfastq_job" "Cellranger mkfastq"
        sleep 10
        # check_results "$check_list" "$base_fastq_path" "dir" "outs" "Cellranger mkfastq completed!" "Cellranger mkfastq did not complete."
        cleanup_snakemake_logs "$(pwd)" "$log_folder" "$base_log_path"
    else
        echo -e "\nRequired Cellranger mkfastq results already exist."
    fi
    chmod -R 777 "$base_fastq_path" >/dev/null 2>&1 
fi


#---------------------------------------------------------------------------------------------------------------------------------------
###### Submit Cellranger counts jobs
if [ "$run_RNAcounts" = "true" ]; then
    # make custom reference data if custom genes are provided
    custom_csv="$base_sheet_folder/custom_refdata.csv"
    if [[ -f "$custom_csv" ]]; then
        while IFS=',' read -r sample idx specie gene_name; do
            bash "$base_src_path/Cellranger/make_refdata.sh" "$gene_name" "$specie" "$REF_PATH" "$PKG_PATH"
        done < "$custom_csv"
    fi
    # make RNAcounts input list
    sample_list=$(make_RNAcount_input_list "$base_count_path" "$log_folder")
    sample_list_check=($(echo "${sample_list[@]}" | grep -v '^\s*$'))
    if [ "${#sample_list_check[@]}" -gt 0 ]; then
        # Merge fastq files if additional BCLs exist
        if [[ -f "$base_sheet_folder/merge.csv" ]]; then

            # Download additional fastq files from google cloud if they do not exist locally
            add_pairs=$(awk -F ',' 'NR>1 {print $2 "," $3; print $4 "," $5}' "$base_sheet_folder/merge.csv" | sort -u)
            if [ -n "$add_pairs" ]; then
                (
                IFS=$'\n'  
                for pair in $add_pairs; do
                    add_bcl=$(echo "$pair" | cut -d ',' -f 1) 
                    add_idx=$(echo "$pair" | cut -d ',' -f 2)
                    add_fastq_path="$BASE_DATA_PATH/$add_bcl/fastq"
                    if [[ ! -z "$add_idx" ]]; then
                        file_exist=false
                        if find "${add_fastq_path}/${add_idx}/" -name "${add_idx}*_S*_L*_R1_*" -print -quit | grep -q . && find "${add_fastq_path}/${add_idx}/" -name "${add_idx}*_S*_L*_R2_*" -print -quit | grep -q .; then
                            file_exist=true
                        fi
                        if [ "$download_fastq" = "true" ] || [ "$file_exist" = "false" ]; then
                            download_files "$add_bcl" "$add_fastq_path" "" "$add_idx" "$base_log_path" "$PKG_PATH" "$GOOGLE_CLOUD_BUCKET" "$use_cluster" "true"
                        fi
                    fi
                done
                )
            fi
            # Link fastq files to the base fastq path
            bash "$base_src_path/Cellranger/merge_fastq.sh" "$base_sheet_folder" "rna"
        fi

        # Run cellranger counts
        echo -e "\n------------------------ Running cellranger counts ------------------------ "
        check_list=()
        for item in "${sample_list_check[@]}"; do
            modified_item=$(echo "$item" | sed -e 's/count_//g' -e 's/vdj_//g' -e 's/multiome_//g')
            check_list+=("$modified_item")
            # download files from google cloud if files do not exist or forced to download
            file_exist=false
            if find "${base_fastq_path}/${modified_item}/" -name "${modified_item}*_S*_L*_R1_*" -print -quit | grep -q . && find "${base_fastq_path}/${modified_item}/" -name "${modified_item}*_S*_L*_R2_*" -print -quit | grep -q .; then
                file_exist=true
            fi
            if [ "$download_fastq" = "true" ] || [ "$file_exist" = "false" ]; then
                download_files "$bcl" "$base_fastq_path" "" "$modified_item" "$base_log_path" "$PKG_PATH" "$GOOGLE_CLOUD_BUCKET" "$use_cluster" "true"
            fi
        done

        submit_batch_jobs "$use_cluster" "$sample_list" "$size_cr" "submit_RNAcounts_job" "Cellranger count"
        sleep 10
        # check_results "$check_list" "$base_count_path" "dir" "outs" "Cellranger counts completed!" "Cellranger counts did not complete."
        cleanup_snakemake_logs "$(pwd)" "$log_folder" "$base_log_path"
    else
        echo -e "\nNo samples to run or required Cellranger count results already exist."
    fi
    chmod -R 777 "$base_count_path" >/dev/null 2>&1 
fi


#---------------------------------------------------------------------------------------------------------------------------------------
###### Submit Cellbender jobs
if [ "$run_cellbender" = "true" ]; then
    sample_list=$(make_cellbener_input_list "$base_count_path" "$log_folder")
    sample_list_check=($(echo "${sample_list[@]}" | grep -v '^\s*$'))
    if [ "${#sample_list_check[@]}" -gt 0 ]; then
        echo -e "\n------------------------ Running Cellbender ------------------------ "
        check_list=()
        for item in "${sample_list[@]}"; do
            modified_item=$(echo "$item" | sed -e 's/count_//g' -e 's/multiome_//g')
            check_list+=("$modified_item")
        done
        submit_batch_jobs "$use_cluster" "$sample_list" "$size_cb" "submit_Cellbender_job" "Cellbender" "none" "$expected_cells" "$total_droplets_included"
        sleep 10
        # check_results "$check_list" "$base_count_path" "file" "cellbender_outs/cellbender_output_filtered.h5" "Cellbender completed!" "Cellbender did not complete."
        cleanup_snakemake_logs "$(pwd)" "$log_folder" "$base_log_path"
    else
        echo -e "\nNo samples to run or required Cellbender results already exist."
    fi
    chmod -R 777 "$base_count_path" >/dev/null 2>&1 
    wait
fi


#---------------------------------------------------------------------------------------------------------------------------------------
###### Submit SBcounts jobs 
if [ "$run_SBcounts" = "true" ]; then
    # Merge fastq files if additional BCLs exist
    if [[ -f "$base_sheet_folder/merge.csv" ]]; then

        # Download additional fastq files from google cloud if they do not exist locally
        add_pairs=$(awk -F ',' 'NR>1 {print $2 "," $3; print $4 "," $5}' "$base_sheet_folder/merge.csv" | sort -u)
        if [ -n "$add_pairs" ]; then
            (
            IFS=$'\n'  
            for pair in $add_pairs; do
                add_bcl=$(echo "$pair" | cut -d ',' -f 1) 
                add_idx=$(echo "$pair" | cut -d ',' -f 2)
                add_fastq_path="$BASE_DATA_PATH/$add_bcl/fastq"

                if [[ ! -z "$add_idx" ]]; then
                    file_exist=false
                    if find "${add_fastq_path}/${add_idx}/" -name "${add_idx}*_S*_L*_R1_*" -print -quit | grep -q . && find "${add_fastq_path}/${add_idx}/" -name "${add_idx}*_S*_L*_R2_*" -print -quit | grep -q .; then
                        file_exist=true
                    fi
                    if [ "$download_fastq" = "true" ] || [ "$file_exist" = "false" ]; then
                        download_files "$add_bcl" "$add_fastq_path" "" "$add_idx" "$base_log_path" "$PKG_PATH" "$GOOGLE_CLOUD_BUCKET" "$use_cluster" "true"
                    fi
                fi
            done
            )
        fi

        # Link fastq files to the base fastq path
        bash "$base_src_path/Cellranger/merge_fastq.sh" "$base_sheet_folder" "sb"
    fi
    # Run SBcounts
    readarray -t output < <(make_SBcount_input_list "$base_spatial_path" "$log_folder")
    sample_list=(${output[0]})
    match_rna_list=(${output[1]})
    sample_list_check=($(echo "${sample_list[@]}" | grep -v '^\s*$'))

    # download files from google cloud if files do not exist or forced to download
    for item in "${sample_list_check[@]}"; do
        modified_item=$(echo "$item" | sed -e 's/SBcount_//g')

        file_exist=false
        if find "${base_fastq_path}/${modified_item}/" -name "${modified_item}*_S*_L*_R1_*" -print -quit | grep -q . && find "${base_fastq_path}/${modified_item}/" -name "${modified_item}*_S*_L*_R2_*" -print -quit | grep -q .; then
            file_exist=true
        fi

        if [ "$download_fastq" = "true" ] || [ "$file_exist" = "false" ]; then
            download_files "$bcl" "$base_fastq_path" "" "$modified_item" "$base_log_path" "$PKG_PATH" "$GOOGLE_CLOUD_BUCKET" "$use_cluster" "fasle"
        fi
    done

    if [ "${#sample_list_check[@]}" -gt 0 ]; then
        echo -e "\n------------------------ Running SBcounts ------------------------ "
        submit_batch_jobs "$use_cluster" "${sample_list[*]}" "$size_sb" "submit_SBcount_job" "SBcount" "${match_rna_list[*]}"
        sleep 10
        # check_results "$match_rna_list" "$base_spatial_path" "file" "SBcounts/SBcounts.h5" "SBcount completed!" "SBcount did not complete."
        cleanup_snakemake_logs "$(pwd)" "$log_folder" "$base_log_path"
    else
        echo -e "\nNo samples to run or required SBcount results already exist."
    fi
    chmod -R 777 "$base_spatial_path" >/dev/null 2>&1 
fi


#---------------------------------------------------------------------------------------------------------------------------------------
###### Submit Spatial analysis jobs
if [ "$run_spatial" = "true" ]; then
    sample_list=$(make_analysis_input_list "$base_spatial_path" "$log_folder")
    sample_list_check=($(echo "${sample_list[@]}" | grep -v '^\s*$'))

    if [ "$download_bam" = "true" ]; then
        for item in "${sample_list_check[@]}"; do
            download_files "$bcl" "" "$base_count_path" "$item" "$base_log_path" "$PKG_PATH" "$GOOGLE_CLOUD_BUCKET" "$use_cluster" "false"
        done
    fi

    if [ "${#sample_list_check[@]}" -gt 0 ]; then
        echo -e "\n------------------------ Running Spatial analysis ------------------------ "
        submit_batch_jobs "$use_cluster" "$sample_list" "$size_sp" "submit_Spatial_job" "Spatial analysis"
        sleep 10
        # check_results "$sample_list" "$base_spatial_path" "file" "Positions/seurat.qs" "Spatial analysis completed!" "Spatial analysis did not complete."
        cleanup_snakemake_logs "$(pwd)" "$log_folder" "$base_log_path"
    else
        echo -e "\nNo samples to run or required Spatial analysis results already exist."
    fi
    chmod -R 777 "$base_spatial_path" >/dev/null 2>&1
fi


#---------------------------------------------------------------------------------------------------------------------------------------
###### Move data to the final folder and clean up
if [ "$run_RNAcounts" = "true" ] || [ "$run_cellbender" = "true" ] || [ "$run_spatial" = "true" ] || [ "$mv_file" = "true" ]; then
    submit_move_data_job "$use_cluster"
    cleanup_snakemake_logs "$(pwd)" "$log_folder" "$base_log_path"
fi