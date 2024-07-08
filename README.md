
# Slide-tags

This is a documentation for the [Slide-tags](https://www.nature.com/articles/s41586-023-06837-4) pipeline (Snakemake-based). It is designed to run on both cluster (UGE 8.5.5) and local environment, and is able to process the data from raw BCL files to spatial analysis. 


## 1. **Installation**

Edit `main_env.yml` and `cellbender_env.yml` files located in `workflow/config` folder. <br>
Change the `name` and `prefix` to the path of the conda environment you want to create. 

```bash
conda env create -f main_env.yml
conda env create -f cellbender_env.yml
```


## 2. **Pipeline Structure**
The pipeline is structured as follows:
- **`env`** <br>
    Folder to store the conda environments installed for the pipeline.

- **`pkgs`** <br>
    Folder to store the packages used in the pipeline, including:
    - bcl2fastq2_v2.20.0
    - CellBender-0.2.2
    - CellBender-0.3.0
    - cellranger-7.2.0
    - cellranger-8.0.1
    - cellranger-arc-2.0.2
    - cellranger-atac-2.1.0
    - google-cloud-sdk

- **`reference`** <br>
    Folder to store the reference genome data used in the pipeline, including:
    - refdata-arc-GRCh38-2020-A      
    - refdata-gex-GRCh38-2024-A
    - refdata-arc-mm10-2020-A    
    - refdata-gex-GRCm39-2024-A
    - refdata-gex-mm10-2020-A
    - refdata-cellranger-vdj-GRCh38  
    - refdata-cellranger-vdj-GRCm38
    - custom_refdata    &nbsp; *(folder for custom reference data)*
    - Index_Info    &nbsp; *(folder for the index information from 10x)*
    - probesets &nbsp; *(folder for the custom probesets)* <br> 

- **`workflow`** <br>
Download the `workflow` from this github repository.

- **`data`** <br>
Each run will create a folder named by `BCL`, which contains the following subfolders:
    - `fastq`<br>
        *Contains FASTQ data, generated by `cellranger mkfastq` or `bcl2fastq`, organized by RNA/SB/ATAC Index.*
    - `count`<br>
        *Contains the count data, generated by `cellranger` or `cellbender`, organized by RNA/ATAC Index.*
      - `outs`
      - `cellbender_outs`
    - `spatial`<br>
        *Contains `SBcounts.h5` and `seurat.qs` files, organized by RNA/SB Index.*
      - `Positions`  
      - `SBcounts`
    - `log`<br>
        *Orginized by the sequence of run, current run_folder recorded in `working`.*
      - `input` *Contains the metadata files for cellranger etc.*
      - `main`  *Contains log files for parsing google sheet and Cluster.*
      - `mkfastq_logs`  *Split by Lane number.*
      - `counts_logs`   *Split by RNA/ATAC Index.*
      - `spatial_logs`  *Split by RNA/SB Index.*


## 3. **Set Configurations**

### 3.1 Change the configurations in `config/config.sh` file to fit your environment.

- **`CLUSTER_PATH`** <br>
    Path to the cluster bin (currently using UGE 8.5.5); set blank if running locally.
- **`CONDA_PATH`** <br>
    Path to the conda bin.
- **`ENV_PATH`** <br>
    Path to the main conda environment for python, julia, and R.
- **`PKG_PATH`** <br>
    Path to the package folder for the pipeline.
- **`BASE_DATA_PATH`** <br>
    Path to the store processed data by pipeline, including outputs of `mkfastq`, `RNAcounts` etc.
- **`BCL_MAIN_PATH`** <br>
    Path to the bcl data; use it as default main path for input BCLs in google sheet.
- **`WORKFLOW_PATH`** <br>
    Path to `workflow` folder, DO NOT contain `workflow` in the path.
- **`GOOGLE_SHEET_ID`** <br>
    Google sheet id for the sample metadata. <br>
    Put your `google_key.json` file in the `workflow/config` folder. <br>
    Here is a [Google sheet demo](https://docs.google.com/spreadsheets/d/1BBsWhvu1bHnhDe-B-3CueJjhJ7JwBmdSVsT9_AWCSTQ/edit?gid=565737114#gid=565737114). 
- **`GOOGLE_CLOUD_BUCKET`** <br>
    Google cloud bucket to store the fastq and bam files.
- **`PUCK_PATH`** <br>
    Path to the puck coordiante csv files .
- **`PUCK_IN`** <br>
    Path to the slide-seq puck barcode files.
- **`REF_PATH`** <br>
    Path to the reference genome data.



###  3.2 Change the parameters of chunk_size in `config/config.sh` file to fit your environment.
Define how many samples to run locally every time:
- `CHUNK_SIZE_MKFASTQ`  = 1
- `CHUNK_SIZE_RNACOUNTS` = 1
- `CHUNK_SIZE_CELLBENDER` = 1
- `CHUNK_SIZE_SBCOUNT` = 10
- `CHUNK_SIZE_POSITION` = 10



###  3.3 Change the parameters of resource usage in `config/smk_config.yaml` file to fit your environment.
```yaml
mem_gb: 128
disk_mb: 8192
runtime_min: 60*60*3
threads: 68
```


## 4. **Run the pipeline**

Export workflow path to the PATH:
```
echo "export PATH="$PATH:/Path/to/slidetag/workflow"" >> ~/.bashrc
source ~/.bashrc
```

Run pipeline after filling sample information in the google sheet:

```
slidetag_pipe.sh -bcl bcl_name -ra
```

Get help message:
```
slidetag_pipe.sh -h
```

```
Usage: Slidetag Pipeline [options]

Required:
  -bcl [value]                  Input BCL name.
Options:
  -h                            Display this help message.
  -ra, -run_all                 Run mkfastq, RNAcounts, SBcounts and Spatial analysis.
  -mk, -run_mkfastq             Run cellranger mkfastq or bcl2fastq.
  -cr, -run_RNAcounts           Run cellranger count or cellranger arc.
  -cb, -run_cellbender          Run cellbender based on Cellranegr count results.
  -sb, -run_SBcounts            Run Spatial beads counts.
  -sp, -run_spatial             Run Spatial analysis for cell positionings.
  -us, -use_sheet               Get input sheets form the current working run.
  -mv, -mv_file                 Move results to store path.
  -gb, -generate_bam            Generate bams when running Cellranger.
  -uf, -upload_fastq            Upload fastqs to google bucket.
  -ub, -upload_bam              Upload bams to google bucket.
  -rf, -rm_fastq                Remove local fastqs.
  -rb, -rm_bam                  Remove local bams.
  -df, -download_fastq          Download fastqs from google bucket.
  -db, -download_bam            Download bams from google bucket.
  -f, -force                    Force to re-run selected jobs.
  -ec [value], -expected_cells  From cellbender parameters.
  -td [value], -total_droplets_included
FALSE or NONE at default for the above parameters.
```
