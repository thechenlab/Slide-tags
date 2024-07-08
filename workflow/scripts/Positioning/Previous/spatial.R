lib_path <- Sys.getenv('R_LIBS', unset = NA)
if (is.na(lib_path)) {
  stop("R_LIBS environment variable is not set.")
} else {
  .libPaths(lib_path)
}
set.seed(42)
library(glue) ; g=glue ; len=length
suppressMessages(library(matrixStats))
suppressMessages(library(ggnewscale))
suppressMessages(library(stringdist))
suppressMessages(library(gridExtra))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(cowplot))
suppressMessages(library(viridis))
suppressMessages(library(ggrastr))
suppressMessages(library(stringr))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(dbscan))
suppressMessages(library(rlist))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(furrr))
suppressMessages(library(rhdf5))
suppressMessages(library(hdf5r))
suppressMessages(library(qpdf))
suppressMessages(library(future))
suppressMessages(library(qs))
options(future.globals.maxSize = 1024 * 1024 * 1024)


r_func_path <- Sys.getenv('R_FUNC', unset = NA)
if (is.na(r_func_path)) {
  stop("R_FUNC environment variable is not set.")
} else {
  source(paste0(r_func_path, '/spatial_functions.R'))
}
base_path <- Sys.getenv('DATA_PATH', unset = NA)
if (is.na(base_path)) {
  stop("DATA_PATH environment variable is not set.")
}
ref_path <- Sys.getenv('REF_PATH', unset = NA)
if (is.na(ref_path)) {
  stop("REF_PATH environment variable is not set.")
}

#### Provide BCL name and optional sample_name in BCL folder
# If only process one sample in BCL foder: Rscript spatial.R BCLname [SAMPLEname]
# If more than one sample nedd to be processed: Rscript spatial.R BCLname [SAMPLEname1, SAMPLEname2]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript spatial.R BCLname [Sample1, Sample2, Sample3]", call. = FALSE)
}
BCLname <- args[1]
SAMPLEnames <- ifelse(length(args) >= 2, args[2], NULL)
ncores <- ifelse(length(args) >= 3, as.numeric(args[3]), 1)

if (!is.null(SAMPLEnames)) {
  if (!grepl("^\\[.*\\]$", SAMPLEnames)) {
    stop("Invalid format for sample names. Please use the format: '[Sample1, Sample2, Sample3]'")
  }
  SAMPLEnames <- gsub("\\[|\\]", "", SAMPLEnames)  
  SAMPLEnames <- strsplit(SAMPLEnames, ",\\s*")[[1]]
}
if (is.null(SAMPLEnames)) {
  message("Processing all samples in BCL: ", BCLname)
} else {
  message("Processing samples ", paste(SAMPLEnames, collapse = ", "), " in BCL: ", BCLname)
}

#### default RNA and Spatial storage path
RNApath = paste0(base_path, "/", BCLname, "/count/")
SBpath = paste0(base_path, "/", BCLname, "/spatial/")
CBdictpath = paste0(ref_path, "/3M-february-2018.txt")


#### load filtered_feature_bc_matrix.h5 or cellbender_output_filtered.h5; load molecule_info.h5
subfolders <- list.files(RNApath, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
names(subfolders) <- sapply(subfolders, basename)
if (is.null(SAMPLEnames) | length(SAMPLEnames) == 0) {SAMPLEnames <- basename(subfolders)}
result <- check_and_get_file_paths(SAMPLEnames, RNApath, SBpath)


#### Seurat processing and spatial mapping
for (n in result$name) {
  message(paste0('\nPROCESSING SAMPLE: ', n))
  RNAh5path = result[result$name == n, ]$rna
  Molpath = result[result$name == n, ]$molecule
  SBh5path = result[result$name == n, ]$spatial

  message('\nSEURAT PROCESSING')
  obj <- seurat_process(RNAh5path)
  if (file.exists(Molpath)) {obj <- add_intronic(obj, Molpath)}
  
  message('\nMAPPING SBCOUNT')
  res <- load_puckdf(SBh5path)
  puckdf <- res[[1]]
  obj@misc %<>% append(res[[2]])
  rm(res) ; invisible(gc())
  res <- map_SBcount(obj, SBh5path, CBdictpath)
  obj <- res$obj
  df <- res$df

  message('\nPOSITION ASSIGNMENT FOR NORMAL DBSCAN')
  dbscan_map_res <- run_positioning(df, obj, ncores = 20L)
  obj <- dbscan_map_res$obj
  data.list <- dbscan_map_res$data.list
  coords <- dbscan_map_res$coords
  dt_path <- paste0(gsub(basename(SBh5path), '', SBh5path), 'normal/data/')
  system(paste0("mkdir -p ", dt_path))
  qsave(df, paste0(dt_path, "sb_mapping_df.qs"))
  qsave(data.list, paste0(dt_path, "data.list.qs"))
  qsave(obj, paste0(dt_path, "seurat.qs"))
  write.csv(coords, paste0(dt_path, "coords.csv"), quote=F, row.names=F)
  
  message('\nSUMMARY PLOTS FOR NORMAL DBSCAN')
  fig_path = paste0(gsub(basename(SBh5path), "", SBh5path), "normal/plots")
  plots_summary(RNAh5path, SBh5path, Molpath, fig_path, obj, df, puckdf, coords, data.list)
}
