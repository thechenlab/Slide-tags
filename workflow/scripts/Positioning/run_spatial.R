# set Lib Path
lib_path <- Sys.getenv('R_LIBS', unset = NA)
if (is.na(lib_path)) {stop("R_LIBS environment variable is not set.")}
.libPaths(lib_path)

# load packages
suppressMessages(library(glue)) ; g=glue ; len=length
suppressMessages(library(gridExtra))
suppressMessages(library(magrittr))
suppressMessages(library(jsonlite))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(Seurat))
suppressMessages(library(rlist))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(qpdf))
suppressMessages(library(qs))
suppressMessages(library(rhdf5))
suppressMessages(library(ggrastr))
suppressMessages(library(Matrix))
suppressMessages(library(dbscan))
suppressMessages(library(rdist))

# Suppress warnings
options(warn = -1)
# # Restore warning settings
# options(warn = getOption("warn"))
set.seed(42)
options(future.globals.maxSize = 1024 * 1024 * 1024)


# Get the path to the R function
r_func_path <- Sys.getenv('R_FUNC', unset = NA)
if (is.na(r_func_path)) {stop("R_FUNC environment variable is not set.")}
source(paste0(r_func_path, '/Functions/run_func.R'))
run_positioning_src <- paste0(r_func_path, '/run_positioning.R')


################ Main ################
args <- commandArgs(trailingOnly = TRUE)
result <- check_args_and_get_paths(args)
for (n in result$name) {
  RNAh5path = result[result$name == n, ]$rna
  Molpath = result[result$name == n, ]$molecule
  SummaryPath = result[result$name == n, ]$summary
  SBh5path = result[result$name == n, ]$spatial
  OutPath = result[result$name == n, ]$output
  ncores = result[result$name == n, ]$ncores
  system(g("Rscript {run_positioning_src} {RNAh5path} {Molpath} {SummaryPath} {SBh5path} {OutPath} {ncores}"))
}
