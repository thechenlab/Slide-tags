# set Lib Path
lib_path <- Sys.getenv('R_LIBS', unset = NA)
if (is.na(lib_path)) {stop("R_LIBS environment variable is not set.")}
.libPaths(lib_path)

# load packages
suppressMessages(library(glue)) ; g=glue ; len=length
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(rlist))
suppressMessages(library(rhdf5))
suppressMessages(library(Matrix))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrastr))
suppressMessages(library(cowplot))
suppressMessages(library(viridis))
suppressMessages(library(magrittr))
suppressMessages(library(jsonlite))
suppressMessages(library(gridExtra))
options(warn = -1)


# Get the path to the R function
r_func_path <- Sys.getenv('R_FUNC', unset = NA)
if (is.na(r_func_path)) {stop("R_FUNC environment variable is not set.")}
source(paste0(r_func_path, '/Functions/load_matrix_func.R'))


# parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  sb_path <- args[[1]]
  cb_path <- args[[2]]
  out_path <- "."
} else if (length(args) == 3) {
  sb_path <- args[[1]]
  cb_path <- args[[2]]
  out_path <- args[[3]]
} else {
  stop("Usage: Rscript load_matrix.R SBcounts_path cb_whitelist_path [output_path]", call. = FALSE)
}
if (!dir.exists(out_path)) { dir.create(out_path, recursive = T) }
stopifnot(system(g("h5ls {sb_path}/matrix"), intern=T) %>% strsplit(split = "\\s+") %>% map_chr(pluck(1)) == c("cb_index", "reads", "sb_index", "umi"))


# load the CB whitelist
cb_whitelist = readLines(cb_path)
stopifnot(class(cb_whitelist) == "character")
stopifnot(!duplicated(cb_whitelist))
stopifnot(len(unique(nchar(cb_whitelist))) == 1)
stopifnot(map_lgl(strsplit(cb_whitelist, ""), ~all(. %in% c("A","C","G","T"))))
message(g("{len(cb_whitelist)} cell barcodes loaded"))

# Load the SB count matrix
df = data.frame(cb_index=f("matrix/cb_index"),
                umi_2bit=f("matrix/umi"),
                sb_index=f("matrix/sb_index"),
                reads=f("matrix/reads"))
message(g("{sum(df$reads)} spatial barcode reads loaded"))

# load metadata
metadata = list(SB_info = list(R1s=f("metadata/R1s"),
                               R2s=f("metadata/R2s"),
                               UMI_downsampling=f("metadata/downsampling"),
                               switch_R1R2=f("metadata/switch") %>% as.logical),
                UP_matching = setNames(f("metadata/UP_matching/count"), f("metadata/UP_matching/type")),
                SB_matching = setNames(f("metadata/SB_matching/count"),f("metadata/SB_matching/type")),
                SB_fuzzy_position = setNames(f("metadata/SB_matching/position_count"), f("metadata/SB_matching/position")),
                bead_info = list(num_lowQ=f("metadata/num_lowQbeads"))
)
metadata$SB_filtering = c(reads_total=f("metadata/num_reads"),
                          reads_noumi=sum(metadata$UP_matching[c("umi_N","umi_homopolymer")]),
                          reads_noup=sum(metadata$UP_matching[c("none","GG")]),
                          reads_nosb=sum(metadata$SB_matching[c("none","HD1ambig")]))



# Fuzzy match and convert cb_index from a cb_list index to a cb_whitelist index
message("Performing fuzzy cell-barcode matching")
df %<>% fuzzy_matching(f("lists/cb_list"), cb_whitelist)
invisible(gc())

# Remove chimeric reads
message("Removing chimeras")
df %<>% remove_chimeras
invisible(gc())

# Plot raw spatial data
message("Creating barcode rank plots")
suppressMessages(suppressWarnings(plot_rankplots(df, f, out_path)))
invisible(gc())

# load the puck information
message("Loading the puck")
puckdf <- load_puckdf(f)

message("Making beadplot")
suppressMessages(suppressWarnings(plot_beadplot(df, puckdf, out_path)))

# remove reads with a filtered spatial barcode
message("Removing low-quality spatial barcodes")
m = df$sb_index %in% puckdf$sb_index
metadata$SB_filtering %<>% c(reads_lowQsb=sum(df$reads[!m]))
df = df[m,]

# remove reads that didn't match a called cell
message("Removing non-whitelist cells")
metadata$SB_filtering %<>% c(reads_uncalled=sum(df$reads[df$cb_index < 0]))
df %<>% filter(cb_index > 0)

# count umis
message("Counting UMIs")
stopifnot(sum(metadata$SB_filtering)+sum(df$reads) == 2 * metadata$SB_filtering["reads_total"])
metadata$SB_filtering %<>% c(reads_final=sum(df$reads))
df %<>% count_umis
metadata$SB_filtering %<>% c(UMIs_final=sum(df$umi))
invisible(gc())

# join tables
stopifnot(df$sb_index %in% puckdf$sb_index)
df %<>% left_join(puckdf, by="sb_index") %>% select(cb_index, x, y, umi) %>% arrange(cb_index, desc(umi))
metadata$puck_info$umi_final = map_int(1:len(metadata$puck_info$puck_name), ~filter(df, x >= metadata$puck_info$puck_boundaries[[.]],
                                                                                    x <= metadata$puck_info$puck_boundaries[[.+1]])$umi %>% sum)

# plot metadata
message("Plotting metadata")
plot_metrics(metadata, out_path)

# write output
message("Writing results")
df %>% setNames(c("cb_index","x_um","y_um","umi")) %>% write.table(file.path(out_path, "matrix.csv.gz") %>% gzfile, sep=",", col.names=T, row.names=F, quote=F)
metadata %>% map(as.list) %>% toJSON(pretty=T) %>% writeLines(file.path(out_path, "spatial_metadata.json"))
