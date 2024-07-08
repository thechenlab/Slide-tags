#!/bin/bash
umask 000

idx=$1
reference=$2
vdj_index=$3
chemistry=$4
version=$5
log_file=$6
fastq_dir=$7
multiome_library=$8
atac_count=$9
dir_path=${10}
generate_bam=${11}
config_file=${12}


# Parse the config file
mem_gb=$(grep -A 4 'RNAcounts:' "$config_file" | grep 'mem_gb' | awk '{print $2}')
threads=$(grep -A 4 'RNAcounts:' "$config_file" | grep 'threads' | awk '{print $2}')


# Load cellranger modules
run_cellranger_vdj() {
    echo "Running cellranger vdj for $idx" | tee -a $log_file
    (time stdbuf -oL -eL cellranger vdj \
        --id=$idx \
        --sample=$idx \
        --fastqs=$fastq_dir \
        --reference=$reference \
        --localcores=$threads \
        --localmem=$mem_gb \
        --chain=$vdj_index &>> $log_file || echo "cellranger vdj failed" >> $log_file) &>> $log_file || true
}

run_cellranger_count() {
    echo "Running cellranger count for $idx" | tee -a $log_file
    echo "generate_bam: $generate_bam" &>> $log_file
    if [[ "$generate_bam" = "false" ]]; then
        if [[ "$version" = "V8" ]]; then
            bam_cmd="--create-bam false"
        else
            bam_cmd="--no-bam"
        fi
    fi
    cmd="time stdbuf -oL -eL cellranger count \
        --id=$idx \
        --sample=$idx \
        --fastqs=$fastq_dir \
        --transcriptome=$reference \
        --jobmode=local \
        --disable-ui  \
        --nosecondary \
        --localcores=$threads \
        --localmem=$mem_gb \
        --chemistry=$chemistry \
        --include-introns=true \
        $bam_cmd"
    (eval $cmd &>> $log_file || echo "cellranger count failed" >> $log_file) &>> $log_file || true
}

run_cellranger_arc_count() {
    echo "Running cellranger-arc count for $idx" | tee -a $log_file
    echo "generate_bam: $generate_bam" &>> $log_file
    if [ "$generate_bam" = "false" ]; then
        bam_cmd="--no-bam"
    fi
    cmd="time stdbuf -oL -eL cellranger-arc count \
        --id=$idx \
        --reference=$reference \
        --libraries=$multiome_library \
        --localcores=$threads \
        --localmem=$mem_gb \
        --jobmode=local \
        $bam_cmd"
    (eval $cmd &>> $log_file || echo "cellranger-arc count failed" >> $log_file) &>> $log_file || true
}

run_cellranger_atac_count() {
    echo "Running cellranger-atac count for $idx" | tee -a $log_file
    (time stdbuf -oL -eL cellranger-atac count \
        --id=$idx \
        --fastqs=$fastq_dir \
        --localcores=$threads \
        --localmem=$mem_gb \
        --reference=$reference &>> $log_file || echo "cellranger-atac count failed" >> $log_file) &>> $log_file || true
}

check_and_retry() {
    local func=$1
    if grep -q "RuntimeError: .* is not a pipestance directory" $log_file; then
        echo "Found RuntimeError related to pipestance directory. Re-executing $func..."
        rm -rf "$dir_path"
        sleep 60
        mkdir -p "$dir_path"
        $func
    fi
}


if [[ -n $vdj_index ]]; then
    run_cellranger_vdj
    sleep 30
    check_and_retry run_cellranger_vdj
elif [[ -f $multiome_library ]]; then
    run_cellranger_arc_count
    sleep 30
    check_and_retry run_cellranger_arc_count
elif [[ $atac_count == "Y" ]]; then
    run_cellranger_atac_count
    sleep 30
    check_and_retry run_cellranger_atac_count
else
    run_cellranger_count
    sleep 30
    check_and_retry run_cellranger_count
fi
