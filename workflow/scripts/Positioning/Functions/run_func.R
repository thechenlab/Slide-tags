#########################################################################################
################################ RUN JOBS FUNCTIONS ##################################
#########################################################################################


###### Qiyu added `check_args_and_get_paths` function to read multiple samples
gex_summary_format <- function(summary_path) {
  df <- read.csv(summary_path)
  column_mapping <- c(
    'Estimated.number.of.cells' = 'Estimated Number of Cells',
    'GEX.Mean.raw.reads.per.cell' = 'Mean Reads per Cell',
    'GEX.Median.genes.per.cell' = 'Median Genes per Cell',
    'GEX.Sequenced.read.pairs' = 'Number of Reads',
    'GEX.Valid.barcodes' = 'Valid Barcodes',
    'GEX.Percent.duplicates' = 'Sequencing Saturation',
    'GEX.Q30.bases.in.barcode' = 'Q30 Bases in Barcode',
    'GEX.Q30.bases.in.read.2' = 'Q30 Bases in RNA Read',
    'GEX.Q30.bases.in.UMI' = 'Q30 Bases in UMI',
    'GEX.Reads.mapped.to.genome' = 'Reads Mapped to Genome',
    'GEX.Reads.mapped.confidently.to.genome' = 'Reads Mapped Confidently to Genome',
    'GEX.Reads.mapped.confidently.to.intergenic.regions' = 'Reads Mapped Confidently to Intergenic Regions',
    'GEX.Reads.mapped.confidently.to.intronic.regions' = 'Reads Mapped Confidently to Intronic Regions',
    'GEX.Reads.mapped.confidently.to.exonic.regions' = 'Reads Mapped Confidently to Exonic Regions',
    'GEX.Reads.mapped.confidently.to.transcriptome' = 'Reads Mapped Confidently to Transcriptome',
    'GEX.Reads.mapped.antisense.to.gene' = 'Reads Mapped Antisense to Gene',
    'GEX.Fraction.of.transcriptomic.reads.in.cells' = 'Fraction Reads in Cells',
    'GEX.Total.genes.detected' = 'Total Genes Detected',
    'GEX.Median.UMI.counts.per.cell' = 'Median UMI Counts per Cell'
  )
  extract_info_columns <- names(column_mapping)
  df_extract <- df %>% select(all_of(extract_info_columns))
  colnames(df_extract) <- column_mapping
  
  percentage_columns <- c(
    'Valid Barcodes',
    'Sequencing Saturation',
    'Q30 Bases in Barcode',
    'Q30 Bases in RNA Read',
    'Q30 Bases in UMI',
    'Reads Mapped to Genome',
    'Reads Mapped Confidently to Genome',
    'Reads Mapped Confidently to Intergenic Regions',
    'Reads Mapped Confidently to Intronic Regions',
    'Reads Mapped Confidently to Exonic Regions',
    'Reads Mapped Confidently to Transcriptome',
    'Reads Mapped Antisense to Gene',
    'Fraction Reads in Cells'
  )
  df_extract <- df_extract %>%
    mutate(across(all_of(percentage_columns), ~ paste0(. * 100, "%")))
  
  # compute Sequencing Saturation
  tr <- df$GEX.Sequenced.read.pairs
  ur <- tr * (1 - df$GEX.Percent.duplicates)
  df_extract$`Sequencing Saturation` <- paste0((1 - (ur / tr)) * 100, "%")
  
  # save name
  save_path <- gsub('summary.csv', 'metrics_summary.csv', summary_path)
  write.csv(df_extract, save_path, row.names = FALSE)
}


check_args_and_get_paths <- function(args) {
  # Get base data path
  if (length(args) < 1) {
    stop("Usage: Rscript run_spatial.R BCLname [Sample1, Sample2]", call. = FALSE)
  }
  base_path <- Sys.getenv('DATA_PATH', unset = NA)
  if (is.na(base_path)) {stop("DATA_PATH environment variable is not set.")}
  
  # Parse args
  BCLname <- args[1]
  SAMPLEnames <- ifelse(length(args) >= 2, args[2], NULL)
  ncores <- ifelse(length(args) >= 3, as.numeric(args[3]), 1)

  if (!is.null(SAMPLEnames)) {
    if (!grepl("^\\[.*\\]$", SAMPLEnames)) {
      stop("Invalid format for sample names. Please use the format: '[Sample1, Sample2]'")
    }
    SAMPLEnames <- gsub("\\[|\\]", "", SAMPLEnames)  
    SAMPLEnames <- strsplit(SAMPLEnames, ",\\s*")[[1]]
  }

  # Get basic RNAcounts and SBcounts paths
  RNApath = paste0(base_path, "/", BCLname, "/count/")
  SBpath = paste0(base_path, "/", BCLname, "/spatial/")

  if (!file.exists(RNApath)) {stop("Input BCLname is invalid!")}
  subfolders <- list.files(RNApath, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
  names(subfolders) <- sapply(subfolders, basename)
  if (is.null(SAMPLEnames) | length(SAMPLEnames) == 0) {SAMPLEnames <- basename(subfolders)}
  
  # Initialize a empty table
  files_df <- data.frame(name = character(), 
                         rna = character(), 
                         molecule = character(), 
                         summary = character(),
                         spatial = character(), 
                         output = character(), 
                         ncores = numeric(),
                         stringsAsFactors = FALSE)

  # Save all sample paths
  for (name in SAMPLEnames) {
    cellbender_file <- paste0(RNApath, name, "/cellbender_outs/cellbender_output_filtered.h5")
    rna_file <- paste0(RNApath, name, "/outs/filtered_feature_bc_matrix.h5")
    mol_file <- paste0(RNApath, name, "/outs/molecule_info.h5")
    summary_file <- paste0(RNApath, name, "/outs/metrics_summary.csv")
    gex_summary_file <- paste0(RNApath, name, "/outs/summary.csv")
    gex_mol_file <- paste0(RNApath, name, "/outs/gex_molecule_info.h5")
    sb_file <- paste0(SBpath, name, "/SBcounts/SBcounts.h5")

    rna_path <- ifelse(file.exists(cellbender_file), cellbender_file, 
                       ifelse(file.exists(rna_file), rna_file, NA))
    mol_path <- ifelse(file.exists(mol_file), mol_file, 
                       ifelse(file.exists(gex_mol_file), gex_mol_file, NA))
    summary_path <- ifelse(file.exists(summary_file), summary_file, 
                          ifelse(file.exists(gex_summary_file), gex_summary_file, NA))
    sb_path <- ifelse(file.exists(sb_file), sb_file, NA)

    missing_files <- c()
    if (is.na(rna_path)) {
      missing_files <- c(missing_files, "filtered_feature_bc_matrix.h5")
    }
    if (is.na(mol_path)) {
      missing_files <- c(missing_files, "molecule_info.h5 or gex_molecule_info.h5")
    }
    if (is.na(summary_path)) {
      missing_files <- c(missing_files, "metrics_summary.csv")
    }
    if (is.na(sb_path)) {
      missing_files <- c(missing_files, "SBcounts.h5")
    }
    if (length(missing_files) > 0) {
      stop(paste("Error: Required files are missing:", paste(missing_files, collapse=", "), "."))
    }

    message(paste("\nRNA:", rna_path, "\nMolecule:", mol_path, "\nSpatial:", sb_path))

    out_path <- gsub("SBcounts", "Positions", gsub(basename(sb_path), '', sb_path))

    if (!file.exists(out_path)) {dir.create(out_path, recursive = TRUE)}
    if (summary_path == gex_summary_file) {
      gex_summary_format(summary_path)
    }
    files_df <- rbind(
      files_df, 
      data.frame(name = name, 
                 rna = rna_path, 
                 molecule = mol_path, 
                 summary = summary_path,
                 spatial = sb_path, 
                 output = out_path,
                 ncores = ncores)
    )
  }
  return(files_df)
}



#########################################################################################
### Qiyu modified `load_seurat`and `ReadCB_h5` function to read Cellbender.h5 files
ReadCB_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = "r")
  genomes <- names(x = infile)
  
  if ("matrix" %in% genomes) {
    genome <- "matrix"
  } else {
    stop("matrix slot not in Cellbender.h5 file")
  }
  if (hdf5r::existsGroup(infile, "matrix")) {
    if (use.names) {
      feature_slot <- "features/name"
    }
    else {
      feature_slot <- "features/id"
    }
  }
  else {
    if (use.names) {
      feature_slot <- "gene_names"
    }
    else {
      feature_slot <- "genes"
    }
  }
  counts <- infile[[paste0(genome, "/data")]]
  indices <- infile[[paste0(genome, "/indices")]]
  indptr <- infile[[paste0(genome, "/indptr")]]
  shp <- infile[[paste0(genome, "/shape")]]
  features <- infile[[paste0(genome, "/", feature_slot)]][]
  barcodes <- infile[[paste0(genome, "/barcodes")]]
  sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[], 
                             x = as.numeric(x = counts[]), dims = shp[], repr = "T")
  if (unique.features) {features <- make.unique(names = features)}
  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- as.sparse(x = sparse.mat)
  if (infile$exists(name = paste0(genome, "/features"))) {
    types <- infile[[paste0(genome, "/features/feature_type")]][]
    types.unique <- unique(x = types)
    if (length(x = types.unique) > 1) {
      message("\nGenome ", genome, " has multiple modalities, returning a list of matrices for this genome")
      sparse.mat <- sapply(X = types.unique, FUN = function(x) {
        return(sparse.mat[which(x = types == x), ])
      }, simplify = FALSE, USE.NAMES = TRUE)
    }
  }
  infile$close_all()
  return(sparse.mat)
}



load_seurat <- function(RNApath) {
  message("\nRunning Seurat processing")
  if(str_ends(RNApath, "filtered_feature_bc_matrix.h5")) {
    obj_data <- Seurat::Read10X_h5(RNApath)
    if (all(c("Gene Expression", "Peaks") %in% names(obj_data))) {
      obj <- CreateSeuratObject(counts = obj_data$`Gene Expression`)
      obj[["Peaks"]] <- CreateAssayObject(counts = obj_data$Peaks)
    } else {
      obj <- CreateSeuratObject(obj_data)
    }
  } else {
    obj <- ReadCB_h5(RNApath) %>% CreateSeuratObject()
  }
  # Add metadata
  obj[["cb"]] <- map_chr(colnames(obj), ~sub("-[0-9]*$", "", .)) %T>% {stopifnot(!duplicated(.))}
  obj[["cb_index"]] <- 1:ncol(obj)
  obj[["logumi"]] <- log10(obj$nCount_RNA+1)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT-|mt-)")
  
  # Add tech metadata
  if (all(c("Peaks", "RNA") %in% names(obj@assays))) {
    Misc(obj, "tech") <- "Multiome"
  } else if (names(obj@assays) == "RNA") {
    Misc(obj, "tech") <- "RNA"
  } else {
    stop("The required assays 'RNA' or 'Peaks' are not present in obj.")
  }
  
  # PCA, Cluster, and UMAP
  obj %<>% Seurat::NormalizeData(verbose = FALSE) %>%
    Seurat::FindVariableFeatures(verbose = FALSE) %>%
    Seurat::ScaleData(verbose = FALSE) %>%
    Seurat::RunPCA(npcs=50, verbose = FALSE) 
  gc()
  
  obj %<>% Seurat::FindNeighbors(dims=1:30, verbose = FALSE) %>%
    Seurat::FindClusters(resolution=round(ncol(obj) / 20000, 2), verbose = FALSE)
  suppressMessages(suppressWarnings(
    obj <- Seurat::RunUMAP(obj, dims=1:30, verbose = FALSE, n.epochs=NULL)
  ))
  
  message(paste0('\nNumber of cells: ', ncol(obj)))
  Misc(obj, "RNApath") <- RNApath

  return(obj)
}



load_intronic <- function(obj, molecule_info_path) {
  library(rhdf5)
  fetch <- function(x){return(h5read(molecule_info_path, x))}
  # Add %intronic
  if ("barcode_idx" %in% h5ls(molecule_info_path)$name) {
    barcodes = fetch("barcodes")
    info = data.frame(barcode=fetch("barcode_idx")+1, umi_type=fetch("umi_type"))
    info %<>% group_by(barcode) %>% summarize(numi=n(), pct.intronic=sum(umi_type==0)/numi)
    obj$pct.intronic = info$pct.intronic[match(obj$cb, barcodes[info$barcode])] * 100
  } else {
    print("WARNING: no intronic information found in the molecule_info_path")
  }
  return(obj)
}



# Helper methods
gdraw <- function(text, s=14) {ggdraw()+draw_label(text, size=s)}
plot.tab <- function(df) {return(plot_grid(tableGrob(df,rows=NULL)))}
add.commas <- function(num){prettyNum(num, big.mark=",")}
make.pdf <- function(plots, name, w, h) {
  if ("gg" %in% class(plots) || class(plots)=="Heatmap") {plots = list(plots)}
  suppressMessages(suppressWarnings({
    pdf(file = name, width = w, height = h)
    lapply(plots, function(x) { print(x) })
    garbage <- dev.off()
  }))
}


# Plot metrics summary
plot_metrics_summary <- function(summary_path) {
  if (nchar(summary_path)==0 || !file.exists(summary_path)) {
    make.pdf(gdraw("No metrics_summary.csv found"), file.path(out_path,"RNAmetrics.pdf"), 7, 8)
    return(c())
  }
  plotdf = read.table(summary_path, header=F, sep=",", comment.char="")
  plotdf %<>% t
  rownames(plotdf) = NULL
  
  plot = plot_grid(ggdraw()+draw_label(""),
                   ggdraw()+draw_label(g("10X Metrics Summary")),
                   plot.tab(plotdf),
                   ggdraw()+draw_label(""),
                   ncol=1, rel_heights=c(0.1,0.1,0.7,0.2))
  make.pdf(plot, file.path(out_path, "RNAmetrics.pdf"), 7, 8)
  
  return(setNames(plotdf[,2], plotdf[,1]))
}



# Plot RNA curves
UvsI <- function(obj, molecule_info_path) {
  if (!file.exists(molecule_info_path) || nchar(molecule_info_path) == 0) {
    plot <- gdraw("No molecule_info.h5 found")
    return(plot)
  }
  if (!"barcode_idx" %in% h5ls(molecule_info_path)$name) {
    plot <- gdraw("Unrecognized molecule_info.h5")
    return(plot)
  }
  
  fetch <- function(x){return(h5read(molecule_info_path,x))}
  barcodes = fetch("barcodes")
  molecule_info = data.frame(barcode=fetch("barcode_idx"),
                             umi_type=fetch("umi_type"),
                             reads=fetch("count"))
  
  # Panel 3: downsampling curve
  tab = table(molecule_info$reads)
  downsampling = map_int(seq(0,1,0.05), function(p){
    map2_int(tab, as.numeric(names(tab)), function(v, k){
      length(unique(floor(sample(0:(k*v-1), round(k*v*p), replace=F)/k)))
    }) %>% sum
  })
  plotdf = data.frame(
    x = seq(0,1,0.05)*sum(molecule_info$reads)/1000/1000,
    y = downsampling/1000/1000
  )
  p0 = ggplot(plotdf, aes(x=x,y=y))+geom_line()+theme_bw()+xlab("Millions of reads")+ylab("Millions of filtered UMIs")+ggtitle("RNA Downsampling curve")
  
  df = molecule_info %>% group_by(barcode) %>% summarize(umi=n(), pct.intronic=sum(umi_type==0)/umi) %>% 
    ungroup %>% arrange(desc(umi)) %>% mutate(logumi=log10(umi))
  
  # Panel 2 and 4: intronic density
  if (!is.null(df$pct.intronic) && !all(df$pct.intronic==0)) {
    ct = 500
    if (any(df$umi>=ct)) {
      p1 = df %>% filter(umi>=ct) %>% ggplot(aes(x = logumi, y = pct.intronic)) + 
        geom_bin2d(bins=100) +
        scale_fill_viridis(trans="log", option="A", name="density") + 
        theme_minimal() +
        labs(title = g("Intronic vs. UMI droplets (>{ct} umi)"), x = "logumi", y = "%intronic") + NoLegend()
      
      max_density_x = density(filter(df,umi>=ct,pct.intronic>1/3)$pct.intronic) %>% {.$x[which.max(.$y)]}
      p2 = df %>% filter(umi>=ct) %>% ggplot(aes(x = pct.intronic)) + geom_density() + 
        theme_minimal() + labs(title = g("Intronic density (>{ct} umi)"), x = "%intronic", y = "Density") + 
        geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
        annotate(geom = 'text', label = round(max_density_x, 2), x = max_density_x+0.01, y = Inf, hjust = 0, vjust = 1, col="red")
    } else {
      p1 = gdraw(g("No cells with {ct}+ UMI"))
      p2 = gdraw(g("No cells with {ct}+ UMI"))
    }
  } else {
    p1 = gdraw(g("No intronic information"))
    p2 = gdraw(g("No intronic information"))
  }
  
  # Panel 1: cell barcode knee plot
  df %<>% mutate(index=1:nrow(df), called=barcodes[barcode+1] %in% obj$cb)
  p3 = ggplot(df,aes(x=index,y=umi,col=called))+geom_line()+theme_bw()+scale_x_log10()+scale_y_log10()+
    ggtitle("Barcode rank plot")+xlab("Cell barcodes")+ylab("UMI counts") +
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"))
  
  plot = plot_grid(p3, p1, p0, p2, ncol=2)
  return(plot)
}


# Plot UMAP + metrics
plot_umaps <- function(obj) {
  mytheme <- function(){theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="top", legend.justification="center", legend.key.width=unit(2, "lines"))}
  umap <- DimPlot(obj,label=T) + ggtitle(g("UMAP")) + mytheme() + NoLegend()  + coord_fixed(ratio=1)
  logumi <- VlnPlot(obj,"logumi",alpha=0) + mytheme() + NoLegend() 
  mt <- FeaturePlot(obj,"percent.mt") + ggtitle("%MT") + mytheme() + coord_fixed(ratio=1) + 
    annotate("text", x = Inf, y = Inf, label = g("Median: {round(median(obj$percent.mt), 2)}%\nMean: {round(mean(obj$percent.mt), 2)}%"), hjust=1, vjust=1, size=2.5, color="black")
  if ("pct.intronic" %in% names(obj@meta.data)) {
    intronic <- FeaturePlot(obj,"pct.intronic") + ggtitle("%Intronic") + mytheme() + coord_fixed(ratio=1) +
      annotate("text", x = Inf, y = Inf, label = g("Median: {round(median(obj$pct.intronic), 2)}%\nMean: {round(mean(obj$pct.intronic), 2)}%"), hjust=1, vjust=1, size=2.5, color="black")
  } else {
    intronic <- gdraw("No intronic information")
  }
  
  plot <- plot_grid(umap, logumi, mt, intronic, ncol=2)
  return(plot)
}


# Create DimPlot
plot_clusters <- function(obj, reduction) {
  npucks = (max(obj$x_um,na.rm=T)-min(obj$x_um,na.rm=T))/(max(obj$y_um,na.rm=T)-min(obj$y_um,na.rm=T))
  nclusters = len(unique(obj$seurat_clusters))
  ncols = round(sqrt(npucks*nclusters/2)/npucks*2) 
  
  m = obj@reductions[[reduction]]@cell.embeddings %>% {!is.na(.[,1]) & !is.na(.[,2])}
  title = g("%placed: {round(sum(m)/len(m)*100,2)} ({sum(m)}/{len(m)}) [{reduction}]")
  p1 = DimPlot(obj, reduction=reduction) + coord_fixed(ratio=1) +
    ggtitle(title) + NoLegend() + xlab("x-position (\u00B5m)") + ylab("y-position (\u00B5m)") + 
    theme(axis.title.x=element_text(size=12), axis.text.x=element_text(size=10)) +
    theme(axis.title.y=element_text(size=12), axis.text.y=element_text(size=10))
  p2 = DimPlot(obj, reduction=reduction, split.by="seurat_clusters", ncol=ncols) + theme_void() + coord_fixed(ratio=1) + NoLegend()
  plot = plot_grid(p1, p2, ncol=1, rel_heights=c(0.4,0.6))
  return(plot)
}


# RNA vs SB metrics
plot_RNAvsSB <- function(obj) {

  coords <- Misc(obj, "coords")
  if ("umi_dbscan" %in% names(coords)) {
    obj$sb_umi <- coords$umi_dbscan %>% tidyr::replace_na(0)
  } else if ("umi" %in% names(coords)) {
    obj$sb_umi <- coords$umi %>% tidyr::replace_na(0)
  } else {
    stop("'umi_dbscan' or 'umi' not found in this Seurat object. Please check the field names.")
  }

  obj$clusters <- Misc(obj,"coords")$clusters %>% tidyr::replace_na(0)
  obj$placed <- !is.na(obj$x_um) & !is.na(obj$y_um)
  
  p1 <- ggplot(obj@meta.data, aes(x=log10(nCount_RNA), y=log10(sb_umi), col=placed)) + 
    geom_point(size=0.2) + theme_bw() + xlab("log10 RNA UMI") + ylab("log10 SB UMI") + ggtitle("SB UMI vs. RNA UMI") + 
    labs(color = "placed") +
    theme(legend.position = c(0.95, 0.05),
          legend.justification = c("right", "bottom"),
          legend.background = element_blank(),
          legend.title=element_text(size=10),
          legend.text=element_text(size=8),
          legend.margin=margin(0,0,0,0,"pt"),
          legend.box.margin=margin(0,0,0,0,"pt"),
          legend.spacing.y = unit(0.1,"lines"),
          legend.key.size = unit(0.5, "lines"))
  
  d = obj@meta.data %>% rowwise %>% mutate(x=min(clusters,5)) %>% ungroup
  p2 <- ggplot(d, aes(x=as.factor(x), y=log10(nCount_RNA))) + geom_violin(scale="count") + 
    scale_x_discrete(breaks=min(d$x):max(d$x), labels=(min(d$x):max(d$x)) %>% {ifelse(.==5, "5+", .)}) +
    xlab("DBSCAN clusters") + ylab("log10 RNA UMI") + ggtitle("RNA UMI vs. DBSCAN cluster") + theme_classic()
  
  d = obj@meta.data %>% group_by(seurat_clusters) %>% summarize(pct.placed=paste0(round(sum(placed)/n()*100,2),"%")) %>% setNames(c("cluster","placed"))
  m = ceiling(nrow(d)/2) ; d1 = d[1:m,] ; d2 = d[(m+1):nrow(d),]
  p3 <- plot_grid(plot.tab(d1), plot.tab(d2), ncol=2)
  
  plot = plot_grid(gdraw("RNA vs. SB metrics"),
                   plot_grid(plot_grid(p1,p2,ncol=1), p3, ncol=2, rel_widths=c(0.5,0.5)),
                   ncol=1, rel_heights=c(0.05,0.95))
  return(plot)
}
