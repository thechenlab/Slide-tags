#!/bin/bash
umask 000

gene_name=$1 
species=$2 # human or mouse
ref_path=$3
pkg_path=$4

set -euo pipefail
default_path="${ref_path}"
cellranger_path="${pkg_path}/cellranger-7.2.0/bin/cellranger" 


# Generate external reference data
cd "${default_path}/custom_refdata/${gene_name}"

# Create the new reference data if it does not exist
if [ -d "genome_${species}_${gene_name}" ]; then
    echo "Reference data genome_${species}_${gene_name} already exists."
    exit 0
else
    if [[ -f "${gene_name}.fasta" ]]; then
        echo "Found ${gene_name}.fasta"
        sed "1s/.*/>${gene_name}/" "${gene_name}.fasta" > "${gene_name}.fa"
    elif [[ -f "${gene_name}.fa" ]]; then
        echo "Found ${gene_name}.fa"
    else
        echo "Gene fasta or fa file not found"
        exit 1
    fi

    echo "Creating ${gene_name}.gtf"
    seq_length=$(cat "${gene_name}.fa" | grep -v "^>" | tr -d "\n" | wc -c)
    echo "Gene length: ${seq_length}"
    echo -e "${gene_name}\tunknown\texon\t1\t${seq_length}\t.\t+\t.\tgene_id \"${gene_name}\"; transcript_id \"${gene_name}\"; gene_name \"${gene_name}\"; gene_biotype \"protein_coding\";" > "${gene_name}.gtf"

    # Add the reference data to the cellranger reference
    msg="Adding ${gene_name} to the ${species} reference data"
    if [ $species == "human" ]; then
        if [ ! -f "genome_${species}_add.fa" ]; then
            echo "$msg"
            cp "$default_path/refdata-gex-GRCh38-2020-A/fasta/genome.fa" "genome_${species}_add.fa"
        else 
            echo "genome_${species}_add.fa already exists"
        fi
        if [ ! -f "genes_${species}_add.gtf" ]; then
            echo "$msg"
            cp "$default_path/refdata-gex-GRCh38-2020-A/genes/genes.gtf" "genes_${species}_add.gtf"
        else
            echo "genes_${species}_add.gtf already exists"
        fi

    elif [ $species == "mouse" ]; then
        if [ ! -f "genome_${species}_add.fa" ]; then
            echo "$msg"
            cp "$default_path/refdata-gex-mm10-2020-A/fasta/genome.fa" "genome_${species}_add.fa"
        else
            echo "genome_${species}_add.fa already exists"
        fi
        if [ ! -f "genes_${species}_add.gtf" ]; then
            echo "$msg"
            cp "$default_path/refdata-gex-mm10-2020-A/genes/genes.gtf" "genes_${species}_add.gtf"
        else
            echo "genes_${species}_add.gtf already exists"
        fi
    else
        echo "Species not supported, please choose human or mouse."
        exit 1
    fi

    # check if the gene is already in the reference
    echo "Checking if ${gene_name} is in the reference data of ${species}..."

    # Check and add to FASTA file
    if grep ">" "genome_${species}_add.fa" | grep -q "${gene_name}"; then
        echo "Gene ${gene_name} already in the reference data."
    else
        echo "Gene ${gene_name} not found in the reference data. Adding now..."
        cat "${gene_name}.fa" >> "genome_${species}_add.fa"
        echo "Gene ${gene_name} added to the reference data."
    fi

    # Check and add to GTF file
    if tail "genes_${species}_add.gtf" | grep -q "${gene_name}"; then
        echo "Gene ${gene_name} already in the GTF reference data."
    else
        echo "Gene ${gene_name} not found in the GTF reference data. Adding now..."
        cat "${gene_name}.gtf" >> "genes_${species}_add.gtf"
        echo "Gene ${gene_name} added to the GTF reference data."
    fi

    echo "Creating reference data with ${gene_name}..."
    echo "${cellranger_path} mkref --genome=genome_${species}_${gene_name} --fasta=genome_${species}_add.fa --genes=genes_${species}_add.gtf"
    $cellranger_path mkref --genome=genome_${species}_${gene_name} --fasta=genome_${species}_add.fa --genes=genes_${species}_add.gtf
    sleep 30
    if [ -d "genome_${species}_${gene_name}" ]; then
        echo "Reference data genome_${species}_${gene_name} created successfully."
    else
        echo "Reference data genome_${species}_${gene_name} not created."
        exit 1
    fi
fi