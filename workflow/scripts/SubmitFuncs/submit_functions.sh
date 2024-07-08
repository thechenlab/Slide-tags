#!/bin/bash
umask 000

### monitor job status with qstat when submitting jobs to the cluster
monitor_job() {
    local job_id=$1
    local sleep_time=${2:-300}
    local start_time=$(date +%s)

    while : ; do
        # Capturing both stdout and stderr from qstat to handle errors
        local qstat_output=$(qstat 2>&1)
        # Check if qstat_output contains the critical error and skip further processing in this iteration if so
        if echo "$qstat_output" | grep -q "getgrgid(1015) failed: Numerical result out of range"; then
            echo "Encountered a critical error with qstat. Will retry..."
            sleep 60
            continue
        fi

        # Checking if the job is found in the qstat output
        local job_info=$(echo "$qstat_output" | grep "$job_id" || echo "Job not found")
        if [[ "$job_info" == "Job not found" ]]; then
            local end_time=$(date +%s)
            local elapsed_time=$((end_time - start_time))
            echo "Job $job_id is not in the queue. Total time elapsed: $((elapsed_time / 60)) minutes."
            break
        fi

        # Extracting job status if job is still in the queue
        local current_status=$(echo "$job_info" | awk '{print $5}')
        local current_time=$(date +%s)
        local time_diff=$((current_time - start_time))
        echo -e "\nCurrent status of job: "
        echo -e "$job_id: \n$current_status."
        echo -e "Time elapsed: $((time_diff / 60)) minutes."

        sleep $sleep_time
    done
}


### parse email address and send email notification for the job status
read_sheet_log="$log_folder/main/read_sheet.log"
emj1=' /(j_j)\\ '
emj2=' \\(^o^)/ '
get_email() {
    file_path="$read_sheet_log"
    local email=""  
    while IFS= read -r line; do
        if [[ "$line" =~ Email:\ (.+) ]]; then
            email="${BASH_REMATCH[1]}"
        fi
    done < "$file_path"
    echo "$email"
}
send_email() {
    local msg="$1"
    local email=$(get_email)

    if [[ ! -z "$email" ]]; then
        echo "Sending email to $email"
        # $python_path "$base_src_path/send_mail.py" "$email" "$msg"
    else
        echo "No email address found."
    fi
}


### clean up snakemake logs
cleanup_snakemake_logs() {
    local dirs=("$@") 
    for dir in "${dirs[@]}"; do
        local snakemake_dir="${dir}/.snakemake/"
        if [ -d "$snakemake_dir" ]; then
            chmod -R 777 "$snakemake_dir" >/dev/null 2>&1
            rm -rf "$snakemake_dir" >/dev/null 2>&1
        fi
    done
}


### submit load google sheet job to the cluster or local computer
submit_read_sheet_job() {
    local use_cluster=${1:-false}
    local read_check_log="$log_folder/main/load_check.log"
    local read_sheet_log="$log_folder/main/read_sheet.log"
    local read_check_err="$log_folder/main/load_check.err"

    local cluster_cmd=""
    if [[ "$use_cluster" == true ]]; then
        cluster_cmd="--cluster \"qsub -N Input -l h_vmem=18G -l h_rt=00:10:00 -pe smp 1 -o $read_check_log -e $read_check_err\""
    fi

    local snakemake_command="snakemake -s \"$base_smk_path/read_check.smk\" \
        --directory \"$log_folder\" \
        --config base_path=\"$BASE_DATA_PATH\" \
        src_path=\"$base_src_path\" \
        py_path=\"$python_path\" \
        bcl=\"$bcl\" bcl_path=\"$BCL_MAIN_PATH\" \
        sheet_id=\"$GOOGLE_SHEET_ID\" \
        ref_path=\"$REF_PATH\" \
        run_folder=\"$log_folder\" \
        run_mkfastq=\"$run_mkfastq\" \
        run_RNAcounts=\"$run_RNAcounts\" \
        run_cellbender=\"$run_cellbender\" \
        run_SBcounts=\"$run_SBcounts\" \
        run_spatial=\"$run_spatial\" \
        --cores all --quiet -j 3 --notemp \
        --latency-wait 10 \
        --forcerun --nolock --keep-going \
        --configfile \"$smk_config\" \
        $cluster_cmd"
    [ -f "$read_check_err" ] && rm -rf "$read_check_err"
    [ -f "$read_sheet_log" ] && rm -rf "$read_sheet_log"
    [ -f "$read_check_log" ] && rm -rf "$read_check_log"
    eval $snakemake_command

    chomod -R 777 "$log_folder/.snakemake" >/dev/null 2>&1

    # Check if the log file exists and no error in the log file
    sleep 10
    if [ ! -d "$log_folder/input" ] || [ ! -f "$read_sheet_log" ] || grep -q "ERROR" "$read_sheet_log"; then
        echo -e "ERROR: Failed to parse inputs. \nPlease check information in google sheet. \nPlease refer to read_sheet.log for detailed error information."
        exit 1
    else
        echo -e "\nGenerate BAM: $generate_bam" >> "$read_sheet_log"
    fi
    if [ ! -f "$read_check_log" ] || grep -q 'ERROR' "$read_check_log"; then
        echo -e "Error in Checking file dependencies. \nPlease check the following files: \nload_check.log in $(basename $log_folder)"
        exit 1
    fi
}


### Submit jobs in batches to the cluster or local computer
submit_batch_jobs() {
    local use_cluster=$1
    local index_list=$2  
    local chunk_size=$3
    local submit_job=$4  
    local job_type=$5  
    local match_rna_list=${6:-none}    
    local expected_cells=${7:-none}   
    local total_droplets_included=${8:-none}

    if [[ "${index_list[*]}" == "${index_list[0]}" ]]; then
        read -ra index_list <<< "$index_list"
    fi

    if [[ "${match_rna_list[*]}" == "${match_rna_list[0]}" ]]; then
        read -ra match_rna_list <<< "$match_rna_list"
    fi

    declare -a current_jobs

    for ((i=0; i<${#index_list[@]}; i+=chunk_size)); do
        local chunk=("${index_list[@]:$i:$chunk_size}")
        echo -e "\nRunning $job_type for ${chunk[@]}"
        local sample_py_list=$(IFS=,; echo "[${chunk[*]}]")

        if [[ "$use_cluster" == true ]]; then
            # Submit the job based on cluster usage
            if [[ "$job_type" == "Cellbender" ]]; then
                local job_id=$($submit_job "$sample_py_list" "$expected_cells" "$total_droplets_included" "$use_cluster")
            elif [[ "$job_type" == "SBcount" ]]; then
                local job_id=$($submit_job "$sample_py_list" "$match_rna_list" "$use_cluster")
            else
                local job_id=$($submit_job "$sample_py_list" "$use_cluster")
            fi
            
            if [ -n "$job_id" ]; then
                echo "Submitted $sample_py_list with iob ID: "
                echo "$job_id"
                current_jobs+=("$job_id")
            else
                echo "Failed to submit job for $sample_py_list."
            fi
        else
            # Submit the job locally
            if [[ "$job_type" == "Cellbender" ]]; then
                $submit_job "$sample_py_list" "$expected_cells" "$total_droplets_included" "$use_cluster"
            elif [[ "$job_type" == "SBcount" ]]; then
                $submit_job "$sample_py_list" "$match_rna_list" "$use_cluster"
            else
                $submit_job "$sample_py_list" "$use_cluster"
            fi
        fi
    done

    # Monitor jobs only if clustering is used
    if [[ "$use_cluster" == true ]]; then
        for job_id in "${current_jobs[@]}"; do
            monitor_job "$job_id" 300
        done
    fi
}


### check the results of the submitted jobs
check_results() {
    local samples=$1
    local base_path=$2
    local check_type=$3 
    local file_or_dir_suffix=$4
    local success_message=$5
    local failure_message=$6

    declare -a finished_list
    declare -a unfinished_list

    for sample in "${samples[@]}"; do
        local target="$base_path/$sample/$file_or_dir_suffix"
        if [[ "$check_type" == "dir" && -d "$target" ]] || [[ "$check_type" == "file" && -f "$target" ]]; then
            finished_list+=("$sample")
        else
            unfinished_list+=("$sample")
        fi
    done

    local error_msg="$failure_message Please check logs in $(basename $log_folder)"
    if [ ${#unfinished_list[@]} -eq 0 ]; then
        echo -e "\n$success_message"
        # send_email "${emj2}\nbcl: ${bcl}\nJob ID: \n${finished_list[@]}\n$success_message"
    else
        echo -e "\n$error_msg"
        # send_email "${emj1}\nbcl: ${bcl}\nJob ID: \n${unfinished_list[@]}\n$error_msg"
        exit 1
    fi
}


### MKFASTQ JOBS FUNCTIONS
# generate the list of mkfastq index for cellranger mkfastq
make_mkfastq_input_list() {
    local base_fastq_path="$1"
    local log_folder="$2"
    local re_run="$force"
    local sheet_path="$log_folder/input/split_mkfastq"
    local prefixes=("ATAC_Index_lane" "Index_lane" "ATAC_lane")
    
    for prefix in "${prefixes[@]}"; do
        if [ ! -d "$sheet_path" ]; then
            echo "ERROR: $sheet_path doesn't exist. Please check read_sheet.log"
            exit 1
        fi
        while IFS= read -r -d $'\0' file; do
            local key=$(basename "$file")
            local number_part=${key#*lane}
            number_part=${number_part%%.*}
            local sample="mkfastq_${number_part}"

            target_path="${base_fastq_path}/${sample}/outs/fastq_path"
            if [ ! -d "$target_path" ]; then
                target_path="${base_fastq_path}/${sample}"
            fi

            # clear all outputs if re-run
            if [[ "$re_run" == true ]]; then
                # rm -rf "${base_fastq_path}/${sample}" 2>/dev/null
                # if [[ $? -ne 0 ]]; then
                #     echo "Failed to delete ${base_fastq_path}/${sample}." >&2
                #     echo "Failed to re-run mkfastq due to permission issue" >&2
                #     exit 0
                # fi
                if ! rm -rf "${base_fastq_path}/${sample}" 2>/dev/null; then
                    echo -e "\nFailed to re-run mkfastq due to permission issue. \nCannot delete ${base_fastq_path}/${sample}" >&2
                    exit 1 
                fi
                sleep 10
            fi

            # make sure all fastqs in mkfastq folder are moved to the correct index folder
            if [ -d "$target_path" ]; then
                local fastq_folder=$(find "$target_path" -maxdepth 1 -type d | grep -vE "(Reports|Stats)" | head -3 | tail -1)
                if [ ! -z "$(ls -A "$fastq_folder")" ]; then
                    local sample_index=$(find "$fastq_folder" -type f -name "*_S*_L*_R*.fastq.gz" -printf "%f\n" | cut -d'_' -f1 | sort -u)
                    for idx in $sample_index; do
                        mkdir -p "$base_fastq_path/$idx"
                        mv "$fastq_folder"/"${idx}"* "$base_fastq_path/$idx/"
                    done
                fi
            fi

            # check if fastqs exist for the all samples
            index_column=($(IFS=$'\n' awk -F',' '{print $NF}' "$file"))
            index_column=("${index_column[@]:1}")

            declare -a index_list
            for idx in "${index_column[@]}"; do
                pattern="${base_fastq_path}/${idx}/${idx}_*_L00${number_part}_*.fastq.gz"
                
                if [[ "$re_run" == true ]]; then
                    # rm -rf "${base_fastq_path}/${idx}" 2>/dev/null
                    # if [[ $? -ne 0 ]]; then
                    #     echo "Failed to delete ${base_fastq_path}/${idx}." >&2
                    #     echo "Failed to re-run mkfastq due to permission issue" >&2
                    #     exit 0
                    # fi
                    if ! rm -rf "${base_fastq_path}/${idx}" 2>/dev/null; then
                        echo -e "\nFailed to re-run mkfastq due to permission issue. \nCannot delete ${base_fastq_path}/${idx}" >&2
                        exit 1 
                    fi
                    sleep 10
                fi

                files_found=$(ls $pattern 2> /dev/null)  
                if [ -z "$files_found" ]; then
                    rm -rf "${base_fastq_path}/${sample}" 2>/dev/null
                    index_list+=("${key%.*}") 
                fi
            done
        done < <(find "$sheet_path" -name "${prefix}*" -print0)
    done
    unique_index_list=($(printf "%s\n" "${index_list[@]}" | sort -u))
    echo "${unique_index_list[@]}"
}


# submit cellranger mkfastq job to the cluster or local computer
submit_mkfastq_job() {
    local sample_py_list=$1
    local use_cluster=${2:-false}
    local mf_id="mf$(date +%d%H%M)"
    local main_err="$log_folder/main/mkfastq_all.err"
    local main_log="$log_folder/main/mkfastq_all.log"

    local cluster_cmd=""
    if [[ "$use_cluster" == true ]]; then
        cluster_cmd="--cluster \"qsub -N $mf_id -r yes -l h_vmem=16G -l h_rt=15:00:00 -binding linear:8 -pe smp 8 -o $main_log -e $main_err\""
    fi

    local snakemake_command="snakemake -s \"$base_smk_path/mkfastq_split.smk\" \
        --directory \"$log_folder\" \
        --config base_path=\"$BASE_DATA_PATH\" pkg_path=\"$PKG_PATH\" \
        bcl=\"$bcl\" run_folder=\"$log_folder\" sample=\"$sample_py_list\" \
        --cores all -j \"$chunk_size\" \
        --immediate-submit --notemp --keep-going \
        --latency-wait 60 --forcerun --quiet --nolock \
        --configfile \"$smk_config\" \
        $cluster_cmd"

    [ -f "$main_log" ] && rm -rf "$main_log"
    [ -f "$main_err" ] && rm -rf "$main_err"
    eval $snakemake_command

    chomod -R 777 "$log_folder/.snakemake" >/dev/null 2>&1

    if [[ "$use_cluster" == true ]]; then
        sleep 60
        local mf_job_id=$(qstat | awk -v mf_id="$mf_id" '$3 ~ mf_id {print $1}')
        if [ -n "$mf_job_id" ]; then
            echo "$mf_job_id"
        else
            echo -e "Failed to submit cellranger mkfastq for $sample_py_list.\n" >&2
            return 1
        fi
    fi
}



### RNAcounts JOBS FUNCTIONS
# generate the list of RNAcounts index for cellranger count
make_RNAcount_input_list() {
    local base_count_path="$1"
    local log_folder="$2"
    local re_run="$force"
    local sheet_path="$log_folder/input/split_counts"
    local prefixes=("count_" "vdj_")

    declare -a sample_list
    for prefix in "${prefixes[@]}"; do
        if [ ! -d "$sheet_path" ]; then
            echo "ERROR: $sheet_path doesn't exist. Please check read_sheet.log"
            exit 1
        fi
        while IFS= read -r -d $'\0' file; do
            local key=$(basename "$file")
            local sample=${key#$prefix} 
            sample=${sample%.*}  
            local target_folder="${base_count_path}/${sample}/outs"

            if [[ "$re_run" == true ]]; then
                # rm -rf "${base_count_path}/${sample}" 2>/dev/null
                # if [[ $? -ne 0 ]]; then
                #     echo "Failed to delete ${base_count_path}/${sample}." >&2
                #     echo "Failed to re-run RNAcounts due to permission issue" >&2
                #     exit 0
                # fi
                if ! rm -rf "${base_count_path}/${sample}" 2>/dev/null; then
                    echo "Failed to re-run RNACounts due to permission issue. \nCannot delete ${base_count_path}/${sample}" >&2
                    exit 1 
                fi
                sleep 10
            fi

            if [ ! -d "$target_folder" ]; then
                sample_list+=("${key%.*}") 
            fi
        done < <(find "$sheet_path" -name "${prefix}*" -print0)
    done
    unique_sample_list=($(printf "%s\n" "${sample_list[@]}" | sort -u))
    echo "${unique_sample_list[@]}"
}

# submit cellranger count job to the cluster or local computer
submit_RNAcounts_job() {
    local sample_py_list=$1
    local use_cluster=${2:-false}
    local cr_id="cr$(date +%d%H%M)"
    local main_err="$log_folder/main/RNAcounts_all.err"
    local main_log="$log_folder/main/RNAcounts_all.log"

    local cluster_cmd=""
    if [[ "$use_cluster" == true ]]; then
        cluster_cmd="--cluster \"qsub -N $cr_id -r yes -l h_vmem=16G -l h_rt=30:00:00 -binding linear:8 -pe smp 8 -o $main_log -e $main_err\""
    fi

    local snakemake_command="snakemake -s \"$base_smk_path/RNAcounts_split.smk\" \
        --directory \"$log_folder\" \
        --config base_path=\"$BASE_DATA_PATH\" \
        pkg_path=\"$PKG_PATH\" \
        ref_path=\"$REF_PATH\" \
        src_path=\"$base_src_path\" \
        bcl=\"$bcl\" run_folder=\"$log_folder\" \
        sample=\"$sample_py_list\" \
        --cores all -j \"$chunk_size\" \
        --immediate-submit --notemp --keep-going \
        --latency-wait 60 --forcerun --quiet --nolock \
        --configfile \"$smk_config\" \
        $cluster_cmd"
    [ -f "$main_log" ] && rm -rf "$main_log"
    [ -f "$main_err" ] && rm -rf "$main_err"
    eval $snakemake_command

    chomod -R 777 "$log_folder/.snakemake" >/dev/null 2>&1

    if [[ "$use_cluster" == true ]]; then
        sleep 30 
        local cr_job_id=$(qstat | awk -v cr_id="$cr_id" '$3 ~ cr_id {print $1}')
        if [ -n "$cr_job_id" ]; then
            echo "$cr_job_id"
        else
            echo -e "Failed to submit cellranger count for $sample_py_list.\n" >&2
            return 1
        fi
    fi
}



### CELLBENDER JOBS FUNCTIONS
# generate the list of cellbender index for cellbender count
make_cellbener_input_list() {
    local base_count_path="$1"
    local log_folder="$2"
    local re_run="$force"
    local sheet_path="$log_folder/input/split_counts"
    local prefixes=("count_")

    declare -a sample_list
    for prefix in "${prefixes[@]}"; do
        if [ ! -d "$sheet_path" ]; then
            echo "ERROR: $sheet_path doesn't exist. Please check read_sheet.log"
            exit 1
        fi
        while IFS= read -r -d $'\0' file; do
            key=$(basename "$file")
            sample=${key#$prefix} 
            sample=${sample%.*}  
            target_file="${base_count_path}/${sample}/cellbender_outs/cellbender_output_filtered.h5"
            
            # clear all outputs if re-run
            if [[ "$re_run" == true ]]; then
                # rm -rf "${base_count_path}/${sample}/cellbender_outs" 2>/dev/null
                # if [[ $? -ne 0 ]]; then
                #     echo "Failed to delete ${base_count_path}/${sample}/cellbender_outs." >&2
                #     echo "Failed to re-run cellbender due to permission issue" >&2
                #     exit 0
                # fi
                if ! rm -rf "${base_count_path}/${sample}/cellbender_outs" 2>/dev/null; then
                    echo -e "\nFailed to re-run Cellbener due to permission issue. \nCannot delete ${base_count_path}/${sample}/cellbender_outs" >&2
                    exit 1 
                fi
            fi

            if [ ! -f "$target_file" ]; then
                sample_list+=("$key")
            fi
        done < <(find "$sheet_path" -name "${prefix}*" -print0)
    done
    unique_sample_list=($(printf "%s\n" "${sample_list[@]}" | sort -u))
    echo "${unique_sample_list[@]}"
}

# submit cellbender job to the cluster or local computer
submit_Cellbender_job() {
    local sample_py_list=$1
    local expected_cells=${2:-none}      
    local total_droplets_included=${3:-none}
    local use_cluster=${4:-false}

    local cb_id="cb$(date +%d%H%M)"
    local main_err="$log_folder/main/cellbender_all.err"
    local main_log="$log_folder/main/cellbender_all.log"

    local cluster_cmd=""
    if [[ "$use_cluster" == true ]]; then
        cluster_cmd="--cluster \"qsub -N $cb_id -r yes -l h_vmem=16G -l h_rt=32:00:00 -binding linear:8 -pe smp 10 -o $main_log -e $main_err\""
    fi

    local snakemake_command="snakemake -s \"$base_smk_path/Cellbender_split.smk\" \
        --directory \"$log_folder\" \
        --config base_path=\"$BASE_DATA_PATH\" \
        env_path=\"$ENV_PATH\" \
        bcl=\"$bcl\" run_folder=\"$log_folder\" \
        sample=\"$sample_py_list\" \
        expected_cells=\"$expected_cells\" \
        total_droplets_included=\"$total_droplets_included\" \
        --cores all -j \"$chunk_size\" \
        --immediate-submit --notemp --keep-going \
        --latency-wait 60 --forcerun --quiet --nolock \
        --configfile \"$smk_config\" \
        $cluster_cmd"
    [ -f "$main_log" ] && rm -rf "$main_log"
    [ -f "$main_err" ] && rm -rf "$main_err"
    eval $snakemake_command

    chomod -R 777 "$log_folder/.snakemake" >/dev/null 2>&1

    if [[ "$use_cluster" == true ]]; then
        sleep 60 
        local cb_job_id=$(qstat | awk -v cb_id="$cb_id" '$3 ~ cb_id {print $1}')
        if [ -n "$cb_job_id" ]; then
            echo "$cb_job_id"
        else
            echo -e "Failed to submit cellbender for $sample_py_list.\n" >&2
            return 1
        fi
    fi
}


# SBcounts JOBS FUNCTIONS
# generate the list of SBcounts index for SBcounts
make_SBcount_input_list() {
    local base_spatial_path="$1"
    local log_folder="$2"
    local re_run="$force"
    local sheet_path="$log_folder/input/split_spatial"
    local prefixes=("SBcount_")

    declare -a sample_list
    declare -a match_rna_list
    for prefix in "${prefixes[@]}"; do
        if [ ! -d "$sheet_path" ]; then
            echo "ERROR: $sheet_path doesn't exist. Please check read_sheet.log"
            exit 1
        fi
        while IFS= read -r -d $'\0' file; do
            key=$(basename "$file")
            sb_index=${key#$prefix}
            sb_index=${sb_index%.*}
            match_tb="$base_sheet_folder/name_to_index.csv"
            sample=$(awk -F, -v sb_idx="$sb_index" '$4 == sb_idx {print $2}' $match_tb)
            target_file="${base_spatial_path}/${sample}/SBcounts/SBcounts.h5"
            puck_file="${base_spatial_path}/${sample}/puck/*.csv"

            # clear all outputs if re-run
            if [[ "$re_run" == true ]]; then
                # rm -rf "$target_file" 2>/dev/null
                # rm -rf "$puck_file" 2>/dev/null
                # if [[ $? -ne 0 ]]; then
                #     echo "Failed to delete $target_file." >&2
                #     echo "Failed to re-run SBcounts due to permission issue" >&2
                #     exit 0
                # fi
                if ! rm -rf "$target_file" "$puck_file" 2>/dev/null; then
                    echo -e "\nFailed to re-run SBcounts due to permission issue. \nCannot delete $target_file" >&2
                    exit 1 
                fi
                sleep 10
            fi

            if [ ! -f "$target_file" ]; then
                sample_list+=("${key%.*}") 
                match_rna_list+=("$sample")
            fi
        done < <(find "$sheet_path" -name "${prefix}*" -print0)
    done
    echo "${sample_list[*]}"
    echo "${match_rna_list[*]}"
}

# submit SBcounts job to the cluster or local computer
submit_SBcount_job() {
    local sample_py_list=$1
    local match_rna_list=$2
    local use_cluster=${3:-false}

    local sb_id="sb$(date +%d%H%M)"
    local main_err="$log_folder/main/SBcounts_all.err"
    local main_log="$log_folder/main/SBcounts_all.log"

    local cluster_cmd=""
    if [[ "$use_cluster" == true ]]; then
        cluster_cmd="--cluster \"qsub -N $sb_id -r yes -l h_vmem=16G -l h_rt=5:00:00 -binding linear:8 -pe smp 8 -o $main_log -e $main_err\""
    fi

    local snakemake_command="snakemake -s \"$base_smk_path/SBcounts_split.smk\" \
        --directory \"$log_folder\" \
        --config base_path=\"$BASE_DATA_PATH\"
        conda_path=\"$CONDA_PATH\" \
        env_path=\"$ENV_PATH\" \
        ref_path=\"$REF_PATH\" \
        src_path=\"$base_src_path\" \
        bcl=\"$bcl\" \
        puck_path=\"$PUCK_PATH\" \
        puck_in=\"$PUCK_IN\" \
        run_folder=\"$log_folder\" \
        sample=\"$sample_py_list\" \
        --cores all -j \"$chunk_size\" \
        --immediate-submit --notemp --keep-going \
        --latency-wait 60 --forcerun \
        --quiet --nolock \
        --configfile \"$smk_config\" \
        $cluster_cmd"
    [ -f "$main_log" ] && rm -rf "$main_log"
    [ -f "$main_err" ] && rm -rf "$main_err"
    eval $snakemake_command

    chomod -R 777 "$log_folder/.snakemake" >/dev/null 2>&1

    if [[ "$use_cluster" == true ]]; then
        sleep 60 
        local sb_job_id=$(qstat | awk -v sb_id="$sb_id" '$3 ~ sb_id {print $1}')
        if [ -n "$sb_job_id" ]; then
            echo "$sb_job_id"
        else
            echo -e "Failed to submit SBcounts for $match_rna_list.\n" >&2
            return 1
        fi
    fi
}



# SPATIAL JOBS FUNCTIONS
# generate the list of Spatial analysis index for Spatial analysis
make_analysis_input_list() {
    local base_spatial_path="$1"
    local log_folder="$2"
    local re_run="$force"
    local sheet_path="$log_folder/input/split_spatial"
    local prefixes=("SBcount_")

    declare -a match_rna_list
    for prefix in "${prefixes[@]}"; do
        if [ ! -d "$sheet_path" ]; then
            echo "ERROR: $sheet_path doesn't exist. Please check read_sheet.log"
            exit 1
        fi
        while IFS= read -r -d $'\0' file; do
            key=$(basename "$file")
            sb_index=${key#$prefix} 
            sb_index=${sb_index%.*} 
            match_tb="$base_sheet_folder/name_to_index.csv"
            sample=$(awk -F, -v sb_idx="$sb_index" '$4 == sb_idx {print $2}' $match_tb)
            # clear all outputs if re-run
            if [[ "$re_run" == true ]]; then
                # rm -rf "${base_spatial_path}/${sample}/Positions" 2>/dev/null
                # if [[ $? -ne 0 ]]; then
                #     echo "Failed to delete ${base_spatial_path}/${sample}/Positions." >&2
                #     echo "Failed to re-run Spatial analysis due to permission issue" >&2
                #     exit 0
                # fi
                if ! rm -rf "${base_spatial_path}/${sample}/Positions" 2>/dev/null; then
                    echo -e "\nFailed to re-run Spatial analysis due to permission issue. \nCannot delete ${base_spatial_path}/${sample}/Positions" >&2
                    exit 1 
                fi
                sleep 10
            fi
            if [ ! -f "${base_spatial_path}/${sample}/Positions/seurat.qs" ]; then
                match_rna_list+=("$sample")
            fi
        done < <(find "$sheet_path" -name "${prefix}*" -print0)
    done
    unique_index_list=($(printf "%s\n" "${match_rna_list[@]}" | sort -u))
    echo "${unique_index_list[@]}"
}

# submit Spatial analysis job to the cluster or local computer
submit_Spatial_job() {
    local match_rna_list=$1
    local use_cluster=${2:-false}

    local sp_id="sp$(date +%d%H%M)"
    local main_err="$log_folder/main/spatial.err"
    local main_log="$log_folder/main/spatial.log"

    local cluster_cmd=""
    if [[ "$use_cluster" == true ]]; then
        cluster_cmd="--cluster \"qsub -N $sp_id -r yes -l h_vmem=12G -l h_rt=5:00:00 -binding linear:8 -pe smp 8 -o $main_log -e $main_err\""
    fi
    
    local snakemake_command="snakemake -s \"$base_smk_path/Spatial_split.smk\" \
        --directory \"$log_folder\" \
        --config base_path=\"$BASE_DATA_PATH\" \
        conda_path=\"$CONDA_PATH\" \
        env_path=\"$ENV_PATH\" \
        ref_path=\"$REF_PATH\" \
        src_path=\"$base_src_path\" \
        bcl=\"$bcl\" \
        rna_list=\"$match_rna_list\" \
        run_folder=\"$log_folder\" \
        --cores all -j \"$chunk_size\" \
        --immediate-submit --notemp --keep-going \
        --latency-wait 60 --forcerun \
        --quiet --nolock \
        --configfile \"$smk_config\" \
        $cluster_cmd"
    [ -f "$main_log" ] && rm -rf "$main_log"
    [ -f "$main_err" ] && rm -rf "$main_err"
    eval $snakemake_command

    chomod -R 777 "$log_folder/.snakemake" >/dev/null 2>&1

    if [[ "$use_cluster" == true ]]; then
        sleep 60 
        local sp_job_id=$(qstat | awk -v sp_id="$sp_id" '$3 ~ sp_id {print $1}')
        if [ -n "$sp_job_id" ]; then
            echo "$sp_job_id"
        else
            echo -e "Failed to submit Spatial analysis for $match_rna_list.\n" >&2
            return 1
        fi
    fi
}


### MOVE DATA FUNCTIONS
submit_move_data_job() {
    local use_cluster=${1:-false}
    local main_log="$log_folder/main/move_data.log"
    local main_err="$log_folder/main/move_data.err"

    echo -e "\n------------------------ Moving results to store path ------------------------ "

    local cluster_cmd=""
    if [[ "$use_cluster" == true ]]; then
        cluster_cmd="--cluster \"qsub -N move_data -l h_vmem=32G -l h_rt=00:20:00 -pe smp 1 -o $main_log -e $main_err\""
    fi
    local snakemake_command="snakemake -s \"$base_smk_path/move_data.smk\" \
        --directory \"$log_folder\" \
        --config base_path=\"$BASE_DATA_PATH\" \
        bcl=\"$bcl\" \
        run_folder=\"$log_folder\" \
        --cores all --quiet -j 3 --notemp --keep-going \
        --latency-wait 30 --forcerun --nolock \
        --configfile \"$smk_config\" \
        $cluster_cmd"
    [ -f "$main_err" ] && rm -rf "$main_err"
    [ -f "$main_log" ] && rm -rf "$main_log"
    eval $snakemake_command

    chomod -R 777 "$log_folder/.snakemake" >/dev/null 2>&1
    sleep 30
}


### UPLOAD DATA FUNCTIONS
upload_fastq() {
    local bcl="$1"
    local base_fastq_path="$2"
    local base_log_path="$3"
    local pkg_path="$4"
    local bucket_name="$5"
    local use_cluster="$6"

    echo -e "\n------------------------ Uploading FASTQ to Google Cloud Bucket ------------------------"

    mkdir -p "$base_log_path/upload"

    if [ "$use_cluster" = "true" ]; then
        local upload_log="$base_log_path/upload/upload_fastq.log"
        local upload_err="$base_log_path/upload/upload_fastq.err"
        [ -f "$upload_log" ] && rm -rf "$upload_log"
        [ -f "$upload_err" ] && rm -rf "$upload_err"
        
        job_id=$(qsub -N "upload" -r yes -l h_vmem=64G -l h_rt=1:00:00 -pe smp 2 -binding linear:2 -o "$upload_log" -e "$upload_err" \
        "$base_src_path/upload_fastq.sh" "$bcl" "$base_fastq_path" "$base_log_path" "$pkg_path" "$bucket_name" | awk '{print $3}')
        echo "Submitted job with ID: $job_id"
        echo "Waiting for uploading to complete..."
        monitor_job "$job_id" 60

        if [ -f "$upload_log" ] && (tail -n 10 "$upload_log" | grep -qE "Finished uploading fastq files.|No fastq folders found"); then
            echo "Uploading FASTQ completed successfully."
        else
            echo "Check the log for errors: $upload_log"
        fi
    else
        bash "$base_src_path/Data2Bucket/upload_fastq.sh" "$bcl" "$base_fastq_path" "$base_log_path" "$pkg_path" "$bucket_name"
    fi
}

upload_bam() {
    local bcl="$1"
    local base_count_path="$2"
    local base_log_path="$3"
    local pkg_path="$4"
    local bucket_name="$5"
    local use_cluster="$6"

    echo -e "\n------------------------ Uploading BAM to Google Cloud Bucket ------------------------"

    mkdir -p "$base_log_path/upload"

    if [ "$use_cluster" = "true" ]; then
        local upload_log="$base_log_path/upload/upload_bam.log"
        local upload_err="$base_log_path/upload/upload_bam.err"
        [ -f "$upload_log" ] && rm -rf "$upload_log"
        [ -f "$upload_err" ] && rm -rf "$upload_err"
        
        job_id=$(qsub -N "upload" -r yes -l h_vmem=64G -l h_rt=1:00:00 -pe smp 2 -binding linear:2 -o "$upload_log" -e "$upload_err" \
        "$base_src_path/Data2Bucket/upload_bam.sh" "$bcl" "$base_count_path" "$base_log_path" "$pkg_path" "$bucket_name" | awk '{print $3}')
        echo "Submitted job with ID: $job_id"
        echo "Waiting for uploading to complete..."
        monitor_job "$job_id" 60

        if [ -f "$upload_log" ] && (tail -n 10 "$upload_log" | grep -qE "Finished uploading bam files.|No bam folders found"); then
            echo "Uploading BAM completed successfully."
        else
            echo "Check the log for errors: $upload_log"
        fi
    else
        bash "$base_src_path/Data2Bucket/upload_bam.sh" "$bcl" "$base_count_path" "$base_log_path" "$pkg_path" "$bucket_name"
    fi
}


### DOWNLOAD DATA FUNCTIONS
download_files() {
    local bcl="$1"
    local base_fastq_path="$2"
    local base_count_path="$3"
    local index_name="$4"
    local base_log_path="$5"
    local pkg_path="$6"
    local bucket_name="$7"
    local use_cluster="$8"
    local index_fastq="$9"

    if [ "$use_cluster" = "true" ]; then
        mkdir -p "$base_log_path/download"
        local download_log="$base_log_path/download/download_files.log"
        local download_err="$base_log_path/download/download_files.err"
        [ -f "$download_log" ] && rm -rf "$download_log"
        [ -f "$download_err" ] && rm -rf "$download_err"
        
        job_id=$(qsub -N "download" -r yes -l h_vmem=64G -l h_rt=0:30:00 -pe smp 2 -binding linear:2 -o "$download_log" -e "$download_err" \
        "$base_src_path/Data2Bucket/download_files.sh" "$bcl" "$base_fastq_path" "$base_count_path" "$index_name" "$pkg_path" "$bucket_name" "$index_fastq" | awk '{print $3}')
        echo "Submitted job with ID: $job_id"
        echo "Waiting for downloading to complete..."
        monitor_job "$job_id" 60

        if [ -f "$download_log" ] && (tail -n 10 "$download_log" | grep -qE "Finished downloading fastq files.|No fastq folders found"); then
            echo "Downloading completed successfully."
        else
            echo "Check the log for errors: $download_log"
        fi
    else
        bash "$base_src_path/Data2Bucket/download_files.sh" "$bcl" "$base_fastq_path" "$base_count_path" "$index_name" "$pkg_path" "$bucket_name" "$index_fastq"
    fi
}