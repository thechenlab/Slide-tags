# set Lib Path
lib_path <- Sys.getenv('R_LIBS', unset = NA)
if (is.na(lib_path)) {stop("R_LIBS environment variable is not set.")}
.libPaths(lib_path)

# load packages
suppressMessages(library(glue)) ; g=glue ; len=length
suppressMessages(library(gridExtra))
suppressMessages(library(magrittr))
suppressMessages(library(Matrix))
suppressMessages(library(jsonlite))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(Seurat))
suppressMessages(library(rlist))
suppressMessages(library(rhdf5))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(qpdf))
suppressMessages(library(qs))
options(warn = -1)


# Get the path to the R function
r_func_path <- Sys.getenv('R_FUNC', unset = NA)
if (is.na(r_func_path)) {stop("R_FUNC environment variable is not set.")}
source(paste0(r_func_path, '/Functions/run_func.R'))
load_matrix_src <- paste0(r_func_path, '/load_matrix.R')
positioning_src <- paste0(r_func_path, '/positioning.R')


# Load arguments
args <- commandArgs(trailingOnly = TRUE)
RNApath <- args[[1]]
molecule_info_path <- args[[2]]
summary_path <- args[[3]]
SBpath <- args[[4]]
out_path <- args[[5]]
ncores <- as.numeric(args[[6]])

### Load the RNA ###############################################################
obj <- load_seurat(RNApath)
if (file.exists(molecule_info_path)) { obj %<>% load_intronic(molecule_info_path) }
Misc(obj, "RNA_metadata") = plot_metrics_summary(summary_path)

### Plots ###############################################################
plot <- UvsI(obj, molecule_info_path)
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path,"RNA.pdf"), 7, 8)))

plot <- plot_umaps(obj)
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path,"UMAP.pdf"), 7, 8)))

### Load the Spatial ###########################################################
# Load the spatial barcode counts matrix and fuzzy match to our called-cells whitelist
message("Generating the matrix...")
cb_whitelist = unname(obj$cb)
writeLines(cb_whitelist, file.path(out_path, "cb_whitelist.txt"))
system(g("Rscript {load_matrix_src} {SBpath} {file.path(out_path, 'cb_whitelist.txt')} {out_path}"))
stopifnot(file.exists(file.path(out_path, "matrix.csv.gz"), file.path(out_path, "spatial_metadata.json")))
Misc(obj, "spatial_metadata") <- fromJSON(file.path(out_path, "spatial_metadata.json"))

# Assign a position to each whitelist cell
message("Positioning cells...")
system(g("Rscript {positioning_src} {file.path(out_path, 'matrix.csv.gz')} {out_path} {ncores}"))
stopifnot(file.exists(file.path(out_path,"coords.csv")))
coords <- read.table(file.path(out_path,"coords.csv"), header=T, sep=",")
coords %<>% right_join(data.frame(cb_index=1:len(cb_whitelist)), by = "cb_index") %>% arrange(cb_index) # fill in cells with no spatial data
Misc(obj, "coords") <- coords

# Create spatial reduction
stopifnot(nrow(coords) == ncol(obj), coords$cb_index == obj$cb_index)
obj$x_um <- coords$x_um
obj$y_um <- coords$y_um

# Add DBSCAN coords
emb = coords %>% select(x_um_dbscan, y_um_dbscan)
colnames(emb) = c("d_1","d_2") ; rownames(emb) = rownames(obj@meta.data)
obj[["dbscan"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "d_", assay = "RNA")

# Add KDE coords
emb = coords %>% mutate(across(everything(), ~ifelse(ratio > 1/3, NA, .))) %>% select(x_um_kde, y_um_kde)
colnames(emb) = c("k_1","k_2") ; rownames(emb) = rownames(obj@meta.data)
obj[["kde"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "k_", assay = "RNA")

# Add KDE-filtered DBSCAN coords
emb = obj@meta.data[,c("x_um","y_um")] ; colnames(emb) = c("s_1","s_2")
obj[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_", assay = "RNA")

# Save the seurat object
qsave(obj, file.path(out_path,"seurat.qs"))

# Plots DBSCAN
plot <- plot_clusters(obj, reduction="dbscan")
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path,"DimPlotDBSCAN.pdf"), 7, 8)))

# Plots KDE
plot <- plot_clusters(obj, reduction="kde")
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path,"DimPlotKDE.pdf"), 7, 8)))

# Plots DimPlot
plot <- plot_clusters(obj, reduction="spatial")
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path,"DimPlot.pdf"), 7, 8)))

# Plots RNA vs SB
plot <- plot_RNAvsSB(obj)
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path, "RNAvsSB.pdf"), 7, 8)))

# Merge the PDF files
plotlist <- c(c("SB.pdf","beadplot.pdf","SBmetrics.pdf"),
              c("DBSCAN.pdf","KDE.pdf","DBSCANvsKDE.pdf","beadplots.pdf"),
              c("RNAmetrics.pdf","RNA.pdf","UMAP.pdf","DimPlot.pdf","DimPlotDBSCAN.pdf","DimPlotKDE.pdf","RNAvsSB.pdf"))
plotorder <- c(8, 9, 10, 1, 2, 4, 5, 6, 11, 12, 13, 14, 3, 7)
suppressMessages({
  pdfs <- file.path(out_path, plotlist[plotorder])
  pdfs %<>% keep(file.exists)
  qpdf::pdf_combine(input=pdfs, output=file.path(out_path,"summary.pdf"))
  file.remove(pdfs)
})
