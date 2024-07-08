# config.sh

# This file contains the configuration for Slide-tag pipeline.
# It should be sourced in the `slidetag_pipe.sh` script.

# # The following variables are defined:
# `CLUSTER_PATH`
#     path to the cluster bin; currently using UGE 8.5.5; set blank if running locally.
# `CONDA_PATH`
#     path to the conda bin.
# `ENV_PATH`
#     path to the main conda environment for python, julia, and R.
# `PKG_PATH`
#     path to the package folder for the pipeline, including:
#         bcl2fastq2_v2.20.0
#         CellBender-0.2.2
#         CellBender-0.3.0
#         cellranger-7.2.0
#         cellranger-8.0.1
#         cellranger-arc-2.0.2
#         cellranger-atac-2.1.0
#         google-cloud-sdk
# `BASE_DATA_PATH`
#     path to the store processed data, including outputs of `mkfastq`, `RNAcounts` etc.
# `BCL_MAIN_PATH`
#     path to the bcl files; use it as default main path for input BCLs in google sheet.
# `WORKFLOW_PATH`
#     path to the `workflow` folder for the pipeline, DO NOT contain `workflow` in the path.
# `GOOGLE_SHEET_ID`
#     the google sheet id for the sample metadata.
# `GOOGLE_CLOUD_BUCKET`
#     the google cloud bucket to store the fastq and bam files.
# `PUCK_PATH`
#     path to the puck coordiante csv files .
# `PUCK_IN`
#     path to the slide-seq puck barcode files.
# `REF_PATH`
#     path to the reference data, including: 
#         refdata-arc-GRCh38-2020-A      
#         refdata-gex-GRCh38-2024-A
#         refdata-arc-mm10-2020-A    
#         refdata-gex-GRCm39-2024-A
#         refdata-gex-mm10-2020-A
#         refdata-cellranger-vdj-GRCh38  
#         refdata-cellranger-vdj-GRCm38
#         custom_refdata                (store custom reference data)
#         Index_Info                    (store the index information from cellranger)
#         probesets                     (store the custom probesets information)

export CLUSTER_PATH=""
export CONDA_PATH=""
export ENV_PATH=""
export PKG_PATH=""
export BASE_DATA_PATH=""
export BCL_MAIN_PATH=""
export WORKFLOW_PATH=""
export GOOGLE_SHEET_ID=""
export GOOGLE_CLOUD_BUCKET=""
export PUCK_PATH=""
export PUCK_IN=""
export REF_PATH=""
export PATH="$CLUSTER_PATH:$PATH"
export PATH="$CONDA_PATH:$PATH"

# The following chunk_size defines the number of samples to run locally every time.
# The default value is 1, which means run one sample at a time.
export CHUNK_SIZE_MKFASTQ=1
export CHUNK_SIZE_RNACOUNTS=1
export CHUNK_SIZE_CELLBENDER=1
export CHUNK_SIZE_SBCOUNT=10
export CHUNK_SIZE_POSITION=10