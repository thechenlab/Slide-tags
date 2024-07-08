
################ load counts info ################
# filtered_feature_bc_matrix.h5 or cellbender_output_filtered.h5; load molecule_info.h5
check_and_get_file_paths <- function(names_list, RNApath, SBpath) {
  if (!file.exists(RNApath)) {message("Input BCLname is invalid!")}
  
  files_df <- data.frame(name=character(), 
                         rna=character(), 
                         molecule=character(), 
                         spatial=character(), 
                         stringsAsFactors=FALSE)
  
  for (name in names_list) {
    cellbender_file <- paste0(RNApath, name, "/cellbender_outs/cellbender_output_filtered.h5")
    rna_file <- paste0(RNApath, name, "/outs/filtered_feature_bc_matrix.h5")
    mol_file <- paste0(RNApath, name, "/outs/molecule_info.h5")
    gex_mol_file <- paste0(RNApath, name, "/outs/gex_molecule_info.h5")
    sb_file <- paste0(SBpath, name, "/SBcounts.h5")
    
    rna_path <- ifelse(file.exists(cellbender_file), cellbender_file, 
                       ifelse(file.exists(rna_file), rna_file, NA))
    mol_path <- ifelse(file.exists(mol_file), mol_file, ifelse(file.exists(gex_mol_file), gex_mol_file, NA))
    sb_path <- ifelse(file.exists(sb_file), sb_file, NA)
    
    missing_files <- c()
    if (is.na(rna_path)) {
      missing_files <- c(missing_files, "filtered_feature_bc_matrix.h5")
    }
    if (is.na(mol_path)) {
      missing_files <- c(missing_files, "molecule_info.h5 or gex_molecule_info.h5")
    }
    if (is.na(sb_path)) {
      missing_files <- c(missing_files, "SBcounts.h5")
    }
    if (length(missing_files) > 0) {
      stop(paste("Error: The following required files are missing:", paste(missing_files, collapse=", "), "."))
    }
    
    message(paste(name, 
                  "\nRNA:", rna_path, 
                  "\nMolecule:", mol_path, 
                  "\nSpatial:", sb_path))
    
    files_df <- rbind(files_df, data.frame(name=name, rna=rna_path, molecule=mol_path, spatial=sb_path))
  }
  
  return(files_df)
}



################ process Seurat ################
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
      message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
      sparse.mat <- sapply(X = types.unique, FUN = function(x) {
        return(sparse.mat[which(x = types == x), ])
      }, simplify = FALSE, USE.NAMES = TRUE)
    }
  }
  infile$close_all()
  return(sparse.mat)
}

seurat_process <- function(RNApath) {
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
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT-|mt-)")
  obj[["logumi"]] <- log10(obj$nCount_RNA+1)
  obj[["cb"]] <- map_chr(colnames(obj), ~sub("-[0-9]*$", "", .))
  
  if (all(c("Peaks", "RNA") %in% names(obj@assays))) {
    Misc(obj, "tech") <- "Multiome"
  } else if (names(obj@assays) == "RNA") {
    Misc(obj, "tech") <- "RNA"
  } else {
    stop("The required assays 'RNA' or 'Peaks' are not present in obj.")
  }

  # PCA, Cluster, and UMAP
  obj %<>% Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(verbose=F) 
  gc()
  
  obj %<>% Seurat::FindNeighbors(dims=1:30) %>%
    Seurat::FindClusters(resolution=round(ncol(obj) / 20000, 2)) %>%
    Seurat::RunUMAP(dims=1:30, verbose=F, n.epochs=NULL)
  
  message(paste0('Number of cells: ', ncol(obj)))
  Misc(obj, "RNApath") <- RNApath
  return(obj)
}

add_intronic<- function(obj, molecule_h5) {
  if (!file.exists(molecule_h5)) {
    message("The molecule_h5 file does not exist. Skipping...")
    return(obj)
  }
  message('add intronic')
  fetch <- function(x){return(h5read(molecule_h5,x))}
  barcodes = fetch("barcodes")
  info = data.frame(barcode=fetch("barcode_idx")+1,
                    feature=fetch("feature_idx")+1,
                    umi=fetch("umi"),
                    umi_type=fetch("umi_type"))
  info %<>% group_by(barcode) %>% summarize(numi=n(), pct.intronic=sum(umi_type==0)/numi)
  obj$pct.intronic = info$pct.intronic[match(obj$cb,barcodes[info$barcode])] * 100
  return(obj)
}



################ process SB matrix ################
# load the puck information
load_puckdf <- function(SBh5path) {
  puckdf = data.frame(sb=h5read(SBh5path, "puck/sb"), 
                      x=h5read(SBh5path, "puck/x"), 
                      y=h5read(SBh5path, "puck/y"), 
                      puck_index=h5read(SBh5path, "puck/puck_index"))
  dups = unique(puckdf$sb[duplicated(puckdf$sb)])
  Ns = puckdf$sb[grepl("N", puckdf$sb)]
  puckdf %<>% filter(!sb %in% dups)
  puckdf %<>% filter(!sb %in% Ns)
  stopifnot(!any(duplicated(puckdf$sb)))
  # puckdfs = map(sort(unique(puckdf$puck_index)), ~filter(puckdf, puck_index==.) %>% select(-puck_index))
  puckdfs <- map(sort(unique(as.character(puckdf$puck_index))), ~filter(puckdf, puck_index == .) %>% select(-puck_index))

  # center the coordinates
  maxs = map_dbl(puckdfs, ~max(.$x))
  mins = map_dbl(puckdfs, ~min(.$x))
  starts = tidyr::replace_na(lag(cumsum(maxs-mins)),0)
  puckdfs %<>% map2(starts, ~mutate(.x,x=x-min(x)+.y)) # line up from y-axis across
  puckdfs %<>% map(~mutate(.,y=y-min(y))) # line up on x-axis
  
  # scale the coordinates (to um)
  num_beads = map_int(puckdfs, nrow)
  scaling_factors = map_dbl(num_beads, get_scaling_factor)
  stopifnot(len(unique(scaling_factors)) == 1)
  scaling_factor = unique(scaling_factors)
  puckdfs %<>% map(~transmute(.,sb=sb, x_um=x*scaling_factor, y_um=y*scaling_factor))
  puckdf = do.call(rbind, puckdfs)
  
  # add sb_index
  sb_list = h5read(SBh5path, "lists/sb_list")
  stopifnot(!any(duplicated(sb_list)))
  puckdf$sb_index = match(puckdf$sb, sb_list)
  puckdf %<>% arrange(sb_index)
  
  puckdf$degen = barcode_degeneracy(puckdf$sb)
  puckdf %<>% select(sb_index, x_um, y_um, degen)
  
  m <- list("puck_name", as.character(h5read(SBh5path, "puck/puck_list")),
            "num_beads", num_beads,
            "xum_start", starts*scaling_factor, 
            "scaling_factor", scaling_factor,
            "dup_beads", length(dups),
            "N_beads", length(Ns),
            "R1s", h5read(SBh5path, "metadata/R1s") %>% map_chr(basename),
            "R2s", h5read(SBh5path, "metadata/R2s") %>% map_chr(basename),
            "switchR1R2", h5read(SBh5path, "metadata/switch") %>% as.logical,
            "UP_matching", setNames(h5read(SBh5path, "metadata/UP_matching/count"), h5read(SBh5path, "metadata/UP_matching/type")),
            "SB_matching", setNames(h5read(SBh5path, "metadata/SB_matching/count"), h5read(SBh5path, "metadata/SB_matching/type")),
            "SB_reads", h5read(SBh5path, "metadata/num_reads"))
  metadata = setNames(m[seq(2, length(m), by=2)], m[seq(1, length(m), by=2)])
  
  return(list(puckdf, metadata))
}
# puck scaling factor
get_scaling_factor <- function(bn) {
  if (bn < 150000) {
    k = 0.73
  } else if (bn < 600000) {
    k = 0.73 * 2
  } else {
    k = 0.645
  }
  return(k)
}
barcode_degeneracy <- function(vec) {
  degenA = str_count(vec, "A")
  degenC = str_count(vec, "C")
  degenG = str_count(vec, "G")
  degenT = str_count(vec, "T")
  degenN = str_count(vec, "N")
  return(pmax.int(degenA, degenC, degenG, degenT, degenN))
}
# remapping
remap_10X_CB <- function(vec) {
  stopifnot(class(vec) == "character")
  basemap = setNames(c("AG","TC","CA","GT"), c("TC","AG","GT","CA"))
  stopifnot(substr(vec,8,9) %in% names(basemap))
  ret = paste0(substr(vec,1,7), basemap[substr(vec,8,9)], substr(vec,10,16))
  stopifnot(length(vec) == length(ret))
  stopifnot(nchar(vec) == nchar(ret))
  return(ret)
}
determine_remap <- function(df, cb_list, cb_whitelist) {
  cbs = cb_list[df$cb_index]
  reads_noremap = df$reads[cbs %in% cb_whitelist] %>% sum
  reads_remap = df$reads[cbs %in% remap_10X_CB(cb_whitelist)] %>% sum
  remap = reads_remap > reads_noremap
  return(remap)
}
listHD1neighbors <- function(input_string) {
  nucleotides <- c('A','C','G','T','N')
  result <- c()
  for (i in 1:nchar(input_string)) {
    for (nuc in nucleotides[nucleotides != substr(input_string, i, i)]) {
      new_string <- paste0(substr(input_string, 1, i-1), nuc, substr(input_string, i+1, nchar(input_string)))
      result <- c(result, new_string)
    }
  }
  return(result)
}
fuzzy_matching <- function(df, cb_list, cb_whitelist) {
  # Check the lists
  stopifnot(!any(duplicated(cb_list)))
  stopifnot(!any(duplicated(cb_whitelist)))
  stopifnot(names(df) == c("cb_index","umi_2bit","sb_index","reads"))
  
  # Remap cb_whitelist
  remap <- determine_remap(df, cb_list, cb_whitelist)
  if (remap) {
    print("Remapping CB whitelist")
    cb_whitelist %<>% remap_10X_CB
  }
  stopifnot(!duplicated(cb_whitelist))
  
  # Exact matching dictionary
  exact_dict = match(cb_list, cb_whitelist)
  
  # HD1 fuzzy matching dictionary
  dict = imap(cb_whitelist, ~data.frame(original=.y, neighbor=listHD1neighbors(.x))) %>% {do.call(rbind,.)}
  dict %<>% filter(neighbor %in% cb_list) # remove barcodes that don't exist
  dict %<>% filter(!neighbor %in% cb_whitelist) # remove exact matches
  HD1ambig = unique(dict$neighbor[duplicated(dict$neighbor)])
  dict %<>% filter(!neighbor %in% HD1ambig) # remove ambiguous matches
  stopifnot(!any(duplicated(dict$neighbor)))
  fuzzy_dict = dict$original[match(cb_list, dict$neighbor)]
  HD1ambig_dict = !is.na(match(cb_list, HD1ambig))
  
  # Perform matching
  df %<>% mutate(exact = exact_dict[cb_index], HD1 = fuzzy_dict[cb_index], HD1ambig = HD1ambig_dict[cb_index])
  stopifnot((!is.na(df$exact)) + (!is.na(df$HD1)) + df$HD1ambig <= 1)
  rm(dict, exact_dict, fuzzy_dict, HD1ambig_dict) ; invisible(gc())
  
  # Record metadata
  CB_matching_type = c("exact", "HD1", "HD1ambig","none")
  CB_matching_count = c(df$reads[!is.na(df$exact)] %>% sum,
                        df$reads[!is.na(df$HD1)] %>% sum,
                        df$reads[df$HD1ambig] %>% sum,
                        df$reads[is.na(df$exact) & is.na(df$HD1) & !df$HD1ambig] %>% sum)
  stopifnot(sum(df$reads) == sum(CB_matching_count))
  
  # Perform the cb_index conversion
  df1 = df %>% filter(!is.na(exact)) %>% mutate(cb_index = exact) %>% select(1:4)
  df2 = df %>% filter(!is.na(HD1)) %>% mutate(cb_index = HD1) %>% select(1:4)
  df3 = df %>% filter(is.na(exact) & is.na(HD1)) %>% mutate(cb_index = -cb_index) %>% select(1:4)
  df2 %<>% group_by(cb_index, umi_2bit, sb_index) %>% summarize(reads=sum(reads)) %>% ungroup
  df12 <- full_join(df1, df2, by = c("cb_index","umi_2bit","sb_index"))
  df12$reads.x %<>% tidyr::replace_na(0) ; df12$reads.y %<>% tidyr::replace_na(0)
  df12 %<>% mutate(reads = reads.x + reads.y) %>% select(-reads.x, -reads.y)
  df = rbind(df12, df3)
  stopifnot(df$cb_index != 0)
  stopifnot(sum(df$reads) == sum(CB_matching_count))
  
  metadata = list(remap=remap, CB_matching=setNames(CB_matching_count, CB_matching_type))
  res = list(df = df, metadata = metadata)
  return(res)
}
# df %<>% group_by(cb_index, sb_index) %>% summarize(umi=n()) 
count_umis <- function(df) {
  stopifnot(names(df) == c("cb_index", "umi_2bit", "sb_index", "reads"))
  gdf = df %>% filter(reads>0) %>% select(cb_index, sb_index) %>% arrange(cb_index, sb_index) 
  bnds = (gdf$cb_index!=lead(gdf$cb_index) | gdf$sb_index!=lead(gdf$sb_index)) %>% tidyr::replace_na(T) %>% which
  gdf %<>% distinct()
  gdf$umi = (bnds - lag(bnds)) %>% tidyr::replace_na(bnds[[1]])
  gdf %<>% arrange(desc(umi))
  return(gdf)
}
# mapping SBcounts
map_SBcount <- function(obj, SBh5path, CBdictpath) {
  # Load the SB counts matrix
  cb_whitelist = unname(obj$cb)
  Misc(obj, "SBh5path") <- SBh5path
  Misc(obj, "called_cells") <- length(cb_whitelist)
  Misc(obj, "RNA_reads") <- sum(obj$nCount_RNA)
  
  if (any(duplicated(cb_whitelist))) {stop("Duplicates found in cb_whitelist.0")}
  if (length(unique(nchar(cb_whitelist))) != 1) {stop("Lengths of barcodes in cb_whitelist are not uniform.")}
  if (!all(map_lgl(strsplit(cb_whitelist, ""), ~all(. %in% c("A", "C", "G", "T"))))) {
    stop("cb_whitelist contains characters other than A, C, G, T.")
  }
  
  df = data.frame(cb_index=h5read(SBh5path, "matrix/cb_index"),
                  umi_2bit=h5read(SBh5path, "matrix/umi"),
                  sb_index=h5read(SBh5path, "matrix/sb_index"),
                  reads=h5read(SBh5path, "matrix/reads"))
  Misc(obj, "SB_reads_filtered") <- sum(df$reads)
  
  # Fuzzy match and convert cb_index from a cb_list index to a cb_whitelist index
  res <- fuzzy_matching(df, h5read(SBh5path, "lists/cb_list"), cb_whitelist)
  df <- res$df ; metadata <- res$metadata
  obj@misc %<>% append(res[[2]])
  rm(res) ; invisible(gc())
  
  # Remove chimeric reads
  message("Removing chimeras")
  df %<>% arrange(cb_index, umi_2bit, desc(reads))
  before_same = tidyr::replace_na(df$cb_index==lag(df$cb_index) & df$umi_2bit==lag(df$umi_2bit), FALSE) 
  after_same = tidyr::replace_na(df$cb_index==lead(df$cb_index) & df$umi_2bit==lead(df$umi_2bit) & df$reads==lead(df$reads), FALSE)
  chimeric = before_same | after_same
  Misc(obj, "SB_reads_filtered_chimeric") <- df[chimeric,]$reads %>% sum
  df = df[!chimeric,] ; rm(chimeric, before_same, after_same)
  
  # remove reads that didn't match a called cell
  df %<>% filter(cb_index > 0)
  
  # Compute metrics
  Misc(obj, "SB_umi_filtered_downsampling") <- h5read(SBh5path, "metadata/downsampling")
  Misc(obj, "SB_umi_final") <- df %>% count_umis %>% pull(umi) %>% sum
  invisible(gc())
  
  return(list(df = df, obj = obj))
}



################ Positioning methods ################
chunk_vector <- function(v, chunk_size) {return(split(v, ceiling(seq_along(v) / chunk_size)))}
# Do a grid search to find the ideal DBSCAN parameters
opt_dbscan <- function(data.list) {
  eps.vec = c(50) ; minPts.vec = c(3:42)
  res = data.frame() ; i = 0
  repeat{
    params = expand.grid(eps.vec,minPts.vec) %>% setNames(c("eps","minPts"))
    row_lists = chunk_vector(1:nrow(params), round(nrow(params)/ncores))
    
    params$pct = furrr::future_map(row_lists, function(v) {
      map_dbl(v, function(i) {
        m = map_lgl(data.list, ~max(dbscan::dbscan(.[c("x_um","y_um")], eps=params$eps[[i]], minPts=params$minPts[[i]], weights=.$umi)$cluster) == 1)
        return(sum(m)/length(m))
      })
    }, .options=furrr_options(seed=T)) %>% flatten_dbl
    
    res = rbind(res, params)
    if (which.max(res$pct)<0.9*nrow(res) || i >= 26) {break}
    minPts.vec = minPts.vec + 40
    i = i + 1
  }
  
  params = res ; rm(res)
  params$is.max = params$pct==max(params$pct)
  eps = params$eps[params$is.max][[1]] ; minPts = params$minPts[params$is.max][[1]]
  pct.placed = round(max(params$pct)*100,2)
  print(g("Optimal eps: {eps}    Optimal minPts: {minPts}    %placed: {pct.placed}"))
  
  return(c(eps,minPts,pct.placed))
}
# Add the DBSCAN clusters to the dataframes
run_dbscan <- function(data.list, eps, minPts) {
  lapply(data.list, function(df){
    df$cluster <- dbscan::dbscan(df[c("x_um","y_um")], eps=eps, minPts=minPts, weights=df$umi)$cluster
    return(df)
  })
}
# assign centroid and record metadata
create_coords <- function(data.list) {
  coords = lapply(data.list, function(df) {
    p = c(x_um=NA,
          y_um=NA,
          DBSCAN_clusters=max(df$cluster),
          SB_umi = sum(df$umi),
          SNR=NA,
          SB_bin = unique(df$bin) %>% {ifelse(is.null(.), NA, .)},
          minPts = unique(df$minPts) %>% {ifelse(is.null(.), NA, .)},
          eps = unique(df$eps) %>% {ifelse(is.null(.), NA, .)},
          pct.placed = unique(df$pct.placed) %>% {ifelse(is.null(.), NA, .)})
    if (max(df$cluster) == 1) {
      sdf = dplyr::filter(df, cluster==1)
      p[["x_um"]] = matrixStats::weightedMedian(sdf$x_um,w=sdf$umi)
      p[["y_um"]] = matrixStats::weightedMedian(sdf$y_um,w=sdf$umi)
      p[["SNR"]] = sum(sdf$umi)/sum(df$umi)
    }
    return(p)
  }) %>% bind_rows %>% as.data.frame %>% mutate(cb_index=as.numeric(names(data.list))) %>% select(cb_index, everything())
  return(coords)
}
# Run these methods on the entire data
normal_positioning <- function(df) {
  stopifnot("umi" %in% colnames(df)) # if this fails, you haven't grouped by cb,sb and counted umis
  data.list = split(df, df$cb_index)
  params = opt_dbscan(data.list)
  data.list %<>% run_dbscan(eps=params[[1]], minPts=params[[2]])
  data.list %<>% map(~mutate(.,eps=params[[1]], minPts=params[[2]], pct.placed=params[[3]]))
  coords <- create_coords(data.list)
  return(list(coords,data.list))
}
# Split the data into 10 SB UMI buckets and run on each
binned_positioning <- function(df) {
  stopifnot("umi" %in% colnames(df)) # if this fails, you haven't grouped by cb,sb and counted umis
  data.list = split(df, df$cb_index)
  
  # create the deciles
  quants = c(0,quantile(map_dbl(data.list,~sum(.$umi)), probs = seq(0.1, 1, by = .1))) %>% unname
  paste("Deciles: ", paste(round(quants),collapse=" "))
  umicounts = map(data.list,~sum(.$umi))
  data.lists=map(1:(len(quants)-1),~data.list[umicounts>quants[.] & umicounts<=quants[.+1]])
  
  # run positioning on each decile
  data.list = map2(data.lists,quants[-1],function(data.list,quant) {
    if (len(data.list) == 0) {print(g({"skipping quantile, empty"}));return(list())}
    params = opt_dbscan(data.list)
    data.list %<>% run_dbscan(eps=params[[1]], minPts=params[[2]])
    data.list %<>% map(~mutate(.,bin=quant,eps=params[[1]],minPts=params[[2]],pct.placed=params[[3]]))
    return(data.list)
  }) %>% list_flatten
  
  stopifnot(len(data.list)==length(unique(df$cb_index)))
  
  coords <- create_coords(data.list)
  return(list(coords, data.list))
}
# Perform positioning at various levels of downsampling
run_positioning <- function(df, obj, ncores = 20L) {
  plan(multisession, workers=ncores)
  original_df <- df
  coords_list = list()
  for (i in seq(1, 1, 0.05)) {
    # Downsample the reads
    # print(g("Downsampling: {round(i*100)}%"))
    df = original_df
    if (i != 1) {
      df %<>% mutate(reads=rmultinom(n=1, size=round(sum(original_df$reads)*i), prob=original_df$reads) %>% as.vector)
    }
    df %<>% count_umis
    # add spatial positions from puck
    df = merge(x=df, y=puckdf, all.x=T, by="sb_index")
    df %<>% filter(!is.na(x_um),!is.na(y_um))
    
    # run normal positioning (see method above)
    res = normal_positioning(df)
    coords = res[[1]] ; data.list = res[[2]] ; rm(res)
    coords_list %<>% list.append(coords)
    gc()
  }
  
  # merge with seurat object
  cb_whitelist = unname(obj$cb)
  coords %<>% mutate(cb = cb_whitelist[cb_index])
  rownames(coords) = paste0(coords$cb,"-1")
  obj = AddMetaData(obj, coords)
  Misc(obj,"pct.placed") = round(sum(!is.na(obj$x_um))/ncol(obj)*100,2)
  emb = obj@meta.data[,c("x_um","y_um")] ; colnames(emb) = c("s_1","s_2")
  obj[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_", assay = "RNA")
  coords %>% select(cb,everything())
  
  return(list(obj = obj, data.list = data.list, coords = coords))
}



################ Create PDFs ################
gdraw <- function(text,s=14) {ggdraw()+ draw_label(text,size=s)}
plot.tab <- function(df) {return(plot_grid(tableGrob(df)))}
add.commas <- function(num){prettyNum(num, big.mark=",")}
make.pdf <- function(plots, name, w, h) {
  if (any(class(plots)=="gg")||class(plots)=="Heatmap") {plots=list(plots)}
  pdf(file=name,width=w,height=h)
  lapply(plots,function(x){print(x)})
  dev.off()
}
UvsI <- function(obj, Molpath) {
  fetch <- function(x){return(h5read(Molpath,x))}
  barcodes = fetch("barcodes")
  molecule_info = data.frame(barcode=fetch("barcode_idx"), umi_type=fetch("umi_type"), reads=fetch("count"))
  
  # Panel 1: downsampling curve
  tab = table(molecule_info$reads)
  downsampling = map_int(seq(0,1,0.05),function(p){
    sum(map2_int(tab, as.numeric(names(tab)), function(v,k){length(unique(floor(sample(0:(k*v-1), round(k*v*p), replace=F)/k)))}))
  })
  plotdf = data.frame(x=seq(0,1,0.05)*sum(molecule_info$reads)/1000/1000, 
                      y=downsampling/1000/1000)
  p0 = ggplot(plotdf, aes(x=x,y=y))+
    geom_line()+theme_bw()+
    xlab("Millions of reads")+
    ylab("Millions of filtered UMIs")+
    ggtitle("RNA Downsampling curve")
  
  df = molecule_info %>% group_by(barcode) %>% 
    summarize(umi=n(), pct.intronic=sum(umi_type==0)/umi) %>% 
    arrange(desc(umi)) %>% 
    mutate(logumi=log10(umi))
  
  # Panel 2 and 3: intronic density
  if (!all(df$pct.intronic==0)) {
    ct = 500
    if (any(df$umi>=ct)) {
      p1 = df %>% filter(umi>=ct) %>% ggplot(aes(x = logumi, y = pct.intronic)) + 
        geom_bin2d(bins=100) +
        scale_fill_viridis(trans="log", option="A", name="density") + 
        theme_minimal() +
        labs(title = g("Intronic vs. UMI droplets (>{ct} umi)"), x = "logumi", y = "%intronic") & NoLegend()
      
      max_density_x = density(filter(df,umi>=ct,pct.intronic>0.35)$pct.intronic) %>% 
        {.$x[which.max(.$y)]}
      p2 = df %>% filter(umi>=ct) %>% 
        ggplot(aes(x = pct.intronic)) +
        geom_density() + 
        theme_minimal() +
        labs(title = g("Intronic density (>{ct} umi)"), x = "%intronic", y = "Density") + 
        geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
        annotate(geom = 'text', label = round(max_density_x, 2), 
                 x = max_density_x+0.01, y = Inf, hjust = 0, vjust = 1, col="red")
    } else {
      p1 = ggdraw()+draw_label(g("No cells with {ct}+ UMI"))
      p2 = ggdraw()+draw_label(g("No cells with {ct}+ UMI"))
    }
  } else {
    p1 = ggdraw()+draw_label("No intronic information")
    p2 = ggdraw()+draw_label("No intronic information")
  }
  
  # Panel 4: cell barcode knee plot
  cb_whitelist = unname(obj$cb)
  df %<>% mutate(index=1:nrow(df), called=barcodes[barcode+1] %in% cb_whitelist)
  p3 = ggplot(df,aes(x=index,y=umi,col=called))+
    geom_line()+theme_bw()+
    scale_x_log10()+
    scale_y_log10()+
    ggtitle("Barcode rank plot")+
    xlab("Cell barcodes")+
    ylab("UMI counts") +
    theme(legend.position = c(0.05, 0.05), 
          legend.justification = c("left", "bottom"), 
          legend.background = element_blank(), 
          legend.spacing.y = unit(0.1,"lines"))
  
  plot = plot_grid(p3,p1,p0,p2,ncol=2)
  return(plot)
}
make_spatial_barcode_plot <- function(df, puckdf, SBh5path) {
  gdf = count_umis(df)
  
  # p1
  cb.data = gdf %>% group_by(cb_index) %>% 
    summarize(umi=sum(umi)) %>% 
    arrange(desc(umi)) %>% 
    {mutate(.,index=1:nrow(.),filter="all cell barcodes")}
  cb.data2 = cb.data %>% filter(cb_index>0) %>% 
    {mutate(.,index=1:nrow(.),filter="called cell barcodes only")}
  sb_pct_in_called_cells = round(sum(filter(cb.data,cb_index>0)$umi)/sum(cb.data$umi)*100,2)
  p1 = ggplot(mapping=aes(x=index, y=umi,col=filter))+
    geom_line(data=cb.data)+geom_line(data=cb.data2) +
    scale_x_log10()+scale_y_log10()+theme_bw()+
    ggtitle("SB UMI per cell")+
    ylab("SB UMI counts")+
    xlab("Cell barcodes") +
    theme(legend.position = c(0.05, 0.05), 
          legend.justification = c("left", "bottom"), 
          legend.background = element_blank(), 
          legend.spacing.y = unit(0.1,"lines"), 
          legend.title=element_blank()) +
    annotate("text", x = Inf, y = Inf, 
             label = g("SB UMI in called cells: {sb_pct_in_called_cells}%"), 
             hjust = 1, vjust = 1.3)
  rm(cb.data, cb.data2) ; invisible(gc())
  
  # p2
  sb.data = gdf %>% group_by(sb_index) %>% 
    summarize(umi=sum(umi)) %>% 
    arrange(desc(umi)) %>% 
    {mutate(.,index=1:nrow(.),filter="all cell barcodes")}
  sb.data2 = gdf %>% filter(cb_index > 0) %>% 
    group_by(sb_index) %>% summarize(umi=sum(umi)) %>% 
    arrange(desc(umi)) %>% 
    {mutate(.,index=1:nrow(.),filter="called cell barcodes only")}
  p2 = ggplot(mapping=aes(x=index,y=umi,col=filter))+
    geom_line(data=sb.data)+
    geom_line(data=sb.data2)+
    scale_x_log10()+
    scale_y_log10()+
    theme_bw()+
    ggtitle("SB UMI per bead")+
    ylab("SB UMI counts")+xlab("Beads")+
    theme(legend.position = c(0.05, 0.05), 
          legend.justification = c("left", "bottom"), 
          legend.background = element_blank(), 
          legend.spacing.y = unit(0.1,"lines"), 
          legend.title=element_blank())
  rm(sb.data, sb.data2) ; invisible(gc())
  
  # p3
  x = seq(0, 1, 0.05) * h5read(SBh5path, "metadata/num_reads")/1000000
  plot.df = data.frame(x=x, y=h5read(SBh5path, "metadata/downsampling")/1000000)
  p3 = ggplot(plot.df, aes(x=x,y=y)) + 
    geom_point() + theme_bw() + 
    xlab("Millions of reads") + 
    ylab("Millions of filtered SB UMIs") + 
    ggtitle("SB downsampling curve")
  
  # p4
  degenmap = setNames(puckdf$degen, puckdf$sb_index)
  gdf.called = filter(gdf, cb_index > 0) %>% 
    group_by(sb_index) %>% 
    summarize(umi=sum(umi)) %>% ungroup %>% 
    mutate(degen = degenmap[as.character(sb_index)]) %>% 
    filter(!is.na(degen)) %>% group_by(degen) %>% 
    summarize(umi=sum(umi)) %>% ungroup %>% 
    arrange(degen) %>% mutate(type="called")
  
  gdf.uncalled = filter(gdf, cb_index < 0) %>% 
    group_by(sb_index) %>% summarize(umi=sum(umi)) %>% 
    ungroup %>% mutate(degen = degenmap[as.character(sb_index)]) %>% 
    filter(!is.na(degen)) %>% group_by(degen) %>% 
    summarize(umi=sum(umi)) %>% ungroup %>% 
    arrange(degen) %>% mutate(type="uncalled")
  plot.df = rbind(gdf.uncalled, gdf.called) %>% 
    mutate(type=factor(type,levels=c("uncalled","called")))
  p4 = ggplot(plot.df, aes(x=degen, y=umi, fill=type)) + 
    geom_col() + theme_bw() +
    scale_x_continuous(breaks = min(plot.df$degen):max(plot.df$degen)) +
    xlab("Spatial barcode degeneracy") + ylab("Number of UMI") + 
    ggtitle("SB degeneracy distribution") +
    theme(legend.position = c(1, 1), legend.justification = c("right", "top"), 
          legend.background = element_blank(), legend.title=element_blank())
  
  plot = plot_grid(p1, p2, p3, p4, ncol=2)
  return(plot)
}
beadplot <- function(sb.data, m, text){
  ggplot(sb.data, aes(x=x_um,y=y_um,col=umi)) +
    rasterize(geom_point(size=0.1), dpi=200) +
    coord_fixed() +
    theme_classic() +
    labs(x="x (\u00B5m)", y="y (\u00B5m)") +
    scale_color_viridis(trans="log", option="B", name="UMI", limits = c(1, m)) + 
    ggtitle(g("SB UMI per bead ({text})"))
}
plot_dbscan <- function(obj, coords) {
  # Panel 1: DBSCAN cluster distribution
  d = data.frame(x=coords$DBSCAN_clusters) %>% 
    rowwise %>% 
    mutate(x=min(x,10)) %>% 
    ungroup
  p1 = ggplot(d,aes(x=x)) + 
    geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))*100), binwidth=.5) +
    geom_text(aes(label = sprintf("%1.0f%%", after_stat(count)/sum(after_stat(count))*100), 
                  y=after_stat(count)/sum(after_stat(count))*100), stat="bin", binwidth=1, vjust=-0.5)+
    theme_classic() + 
    xlab("Num DBSCAN clusters") + 
    ylab("Percent") +
    scale_y_continuous(limits=c(0,100)) +
    scale_x_continuous(breaks=min(d$x):max(d$x), labels=(min(d$x):max(d$x)) %>% {ifelse(.==10, "10+", .)}) +
    ggtitle("DBSCAN cluster distribution")
  
  # Panel 2: SNR density
  max_density_x = density(obj$SNR %>% na.omit) %>% {.$x[which.max(.$y)]}
  max_density_x = median(obj$SNR, na.rm = T)
  p2 = obj@meta.data %>% filter(!is.na(x_um)) %>% ggplot(aes(x = SNR)) +
    geom_density() + 
    theme_minimal() +
    labs(title = "SNR per cell (density)", x = "SNR", y = "Density") + 
    geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
    annotate(geom = 'text', label = round(max_density_x, 2), 
             x = max_density_x+0.01, y = Inf, hjust = 0, vjust = 1, col="red")
  
  # Panel 3: RNA umi vs SB umi
  p3 = data.frame(x=obj$nCount_RNA,y=obj$SB_umi,placed=!is.na(obj$x_um)) %>% 
    {ggplot(.,aes(x=log10(x),y=log10(y),col=placed))+
        geom_point(size=0.2)+theme_bw()+xlab("RNA UMI")+
        ylab("SB UMI")+ggtitle("SB UMI vs. RNA UMI")+
        theme(legend.position = c(0.95, 0.05), 
              legend.justification = c("right", "bottom"), 
              legend.background = element_blank(), 
              legend.title=element_text(size=10), 
              legend.text=element_text(size=8), 
              legend.margin=margin(0,0,0,0,"pt"), 
              legend.box.margin=margin(0,0,0,0,"pt"), 
              legend.spacing.y = unit(0.1,"lines"), 
              legend.key.size = unit(0.5, "lines"))}
  
  # Panel 4: DBSCAN parameters
  df = coords %>% select(SB_bin,minPts,eps,pct.placed) %>% 
    distinct %>% arrange(SB_bin) %>% 
    mutate(SB_bin=round(SB_bin,2), pct.placed=round(pct.placed,2) %>% paste0("%"))
  rownames(df) <- NULL
  p4 = plot_grid(gdraw("DBSCAN parameters"),plot.tab(df),ncol=1,rel_heights=c(1,17))
  
  plot = plot_grid(p1,p2,p3,p4,ncol=2)
  return(plot)
}
metrics_plots <- function(obj) {
  plot.df = list(
    c("Reads",Misc(obj,"SB_reads") %>% add.commas),
    c("Puck file",paste0(Misc(obj,"puck_name") %>% str_replace("Puck_","")%>% str_replace(".csv",""),collapse=", ")),
    c("Number of beads",Misc(obj,"num_beads") %>% add.commas %>% paste0(collapse=", ")),
    c("Scaling factor",Misc(obj,"scaling_factor")),
    c("R1<->R2",Misc(obj,"switchR1R2")),
    c("Remap CB",Misc(obj,"remapCB"))
  ) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("metric","value"))
  p1 = plot_grid(ggdraw()+draw_label("SB metrics"), plot.tab(plot.df), ncol=1, rel_heights=c(0.1,0.9))
  
  UP_matching <- Misc(obj,"UP_matching")
  SB_matching <- Misc(obj,"SB_matching")
  CB_matching <- Misc(obj,"CB_matching")

  plot.df = data.frame(a=c("exact","fuzzy", "none", "GG"),
                       b=c(UP_matching[["-"]],
                           sum(UP_matching[c("1D-","1D-1X","-1X","-1D","-2X")]),
                           UP_matching[["none"]],
                           UP_matching[["GG"]]) %>% {./sum(.)*100} %>% round(2) %>% paste0("%")) %>% arrange(desc(b)) %>% unname
  p2 = plot_grid(ggdraw()+draw_label("UP matching"), plot.tab(plot.df),ncol=1,rel_heights=c(0.1,0.9))
  
  plot.df = data.frame(a=c("exact","fuzzy","none","ambig"),b=SB_matching[c("exact","HD1","none","HD1ambig")] %>% 
                         {./sum(.)*100} %>% round(2) %>% unname %>% paste0("%")) %>% unname
  p3 = plot_grid(ggdraw()+draw_label("SB matching"),plot.tab(plot.df),ncol=1,rel_heights=c(0.1,0.9))
  
  plot.df = data.frame(a=c("exact","fuzzy","none","ambig"),b=CB_matching[c("exact","HD1","none","HD1ambig")] %>% 
                         {./sum(.)*100} %>% round(2) %>% unname %>% paste0("%")) %>% unname
  p4 = plot_grid(ggdraw()+draw_label("CB matching"),plot.tab(plot.df),ncol=1,rel_heights=c(0.1,0.9))
  
  plot.df = list(
    c("Valid UMI",UP_matching[["R1lowQ"]]),
    c("Valid UP",UP_matching[c("none","GG")] %>% sum),
    c("Valid SB",SB_matching[c("none","HD1ambig")] %>% sum),
    c("Valid CB",sum(CB_matching[c("none","HD1ambig")])),
    c("Chimeric",Misc(obj,"SB_reads_filtered_chimeric"))
  ) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("filter","percent"))
  plot.df$percent = as.numeric(plot.df$percent) / (Misc(obj,"SB_reads") - c(0,head(cumsum(plot.df$percent),-1)))
  plot.df[1:4,2] = 1-plot.df[1:4,2] # convert from fraction removed to fraction retained
  plot.df$percent = round(plot.df$percent * 100,2) %>% paste0("%")
  p5 = plot_grid(ggdraw()+draw_label("SB filtering"),plot.tab(plot.df),ncol=1,rel_heights=c(0.1,0.9))
  
  SB_filtering = setNames(plot.df$percent, plot.df$filter)
  
  # Script input parameters (with common prefix removed)
  p6 = ggplot()
  tryCatch( {
    a=Misc(obj,"RNA_path") ; b=Misc(obj,"SB_path")
    m=min(which(!map_lgl(1:min(nchar(a),nchar(b)), ~str_sub(a,1,.)==str_sub(b,1,.))))
    a%<>%str_sub(m-1,999) ; b%<>%str_sub(m-1,999)
    plot.df=data.frame(a=c("RNA_path","SB_path"),b=c(a,b)) %>% unname
  }, error = function(e) {p6 <<- ggdraw()+draw_label("Error") })
  
  plot = plot_grid(
    gdraw("Additional metadata",18),
    plot_grid(p1,p5,ncol=2),
    plot_grid(p2,p3,p4,ncol=3),
    ggdraw()+draw_label(""), #spacer
    ncol=1,
    rel_heights = c(0.27,0.5,0.35,0.25,0.27)
  )
  
  return(list(plot,SB_filtering))
}
sample_bead_plots <- function(data.list, puckdf) {
  plot.sb <- function(subdf) {
    subdf1 <- filter(subdf,cluster==1)
    subdf2 <- filter(subdf,cluster==2)
    ggplot()+coord_fixed(ratio=1,xlim=range(puckdf$x_um),ylim=range(puckdf$y_um))+theme_void()+
      geom_point(data=subdf, mapping=aes(x=x_um,y=y_um,col=umi),size=2,shape=16)+
      geom_point(aes(x=matrixStats::weightedMedian(subdf1$x_um,w=subdf1$umi),
                     y=matrixStats::weightedMedian(subdf1$y_um,w=subdf1$umi)),
                 color="red",shape=0,size=3) + 
      geom_point(aes(x=matrixStats::weightedMedian(subdf2$x_um,w=subdf2$umi),
                     y=matrixStats::weightedMedian(subdf2$y_um,w=subdf2$umi)),
                 color="green",shape=0,size=3) + 
      theme(legend.key.width=unit(0.5,"lines"), 
            legend.position="right", 
            legend.key.height=unit(1,"lines"), 
            legend.title=element_blank(), 
            legend.spacing.y=unit(0.2,"lines"), 
            legend.margin=margin(0,0,0,0,"lines"), 
            legend.box.margin=margin(0,0,0,0,"pt"), 
            legend.box.background=element_blank(), 
            legend.background=element_blank(), 
            legend.direction="vertical", 
            legend.justification="left",
            legend.box.just="left",
            legend.box.spacing=unit(0,"cm"))
  }
  
  list0 = data.list %>% keep(~max(.$cluster)==0) %>% map(~arrange(.,umi) %>% filter(!is.na(x_um), !is.na(y_um)))
  list1 = data.list %>% keep(~max(.$cluster)==1) %>% map(~arrange(.,umi) %>% filter(!is.na(x_um), !is.na(y_um)))
  list2 = data.list %>% keep(~max(.$cluster)==2) %>% map(~arrange(.,umi) %>% filter(!is.na(x_um), !is.na(y_um)))
  
  list0 = list0[as.integer(map(list0,nrow))!=0]
  list1 = list1[as.integer(map(list1,nrow))!=0]
  list2 = list2[as.integer(map(list2,nrow))!=0]
  
  p1 = plot_grid(ggdraw()+draw_label("DBSCAN=0"), map(sample(list0,min(12,len(list0)),replace=F),plot.sb) %>% 
                   {plot_grid(plotlist=.,ncol=3)}, ncol=1,rel_heights=c(0.1,2))
  p2 = plot_grid(ggdraw()+draw_label("DBSCAN=1"), map(sample(list1,min(12,len(list1)),replace=F),plot.sb) %>% 
                   {plot_grid(plotlist=.,ncol=3)}, ncol=1,rel_heights=c(0.1,2))
  p3 = plot_grid(ggdraw()+draw_label("DBSCAN=2"), map(sample(list2,min(12,len(list2)),replace=F),plot.sb) %>% 
                   {plot_grid(plotlist=.,ncol=3)}, ncol=1,rel_heights=c(0.1,2))
  
  plots = list(p1,p2,p3)
  
  return(plots)
}
# plot all figs
plots_summary <- function(RNAh5path, SBh5path, Molpath, fig_path, 
                          obj, df, puckdf, coords, data.list){
  message("Making summary.pdf")
  system(paste0("mkdir ", fig_path))
  
  ###### Page 0: cellranger output
  message('0.cellranger output')
  tryCatch({
    metrics_summary = gsub(basename(RNAh5path), 'metrics_summary.csv', RNAh5path)
    if (!file.exists(metrics_summary)) {
      RNAh5path = gsub("/cellbender_outs/cellbender_output_filtered.h5", "/outs/", RNAh5path)
      metrics_summary = paste0(RNAh5path, 'metrics_summary.csv')
    }
    
    if (file.exists(metrics_summary)) {
      plotdf = read.table(metrics_summary, header=F, comment.char="", sep=",")
      if (nrow(plotdf) == 2) { # count
        plotdf %<>% t()
      } else if (ncol(plotdf) == 6) { # multi
        colnames(plotdf) = as.character(plotdf[1,])
        plotdf = plotdf[-1, c(5,6)]
      }
      rownames(plotdf) = NULL

      # Add metadata to seurat object
      Misc(obj,"RNA_metrics") <- list(plotdf[,1], plotdf[,2])
      
      plot = plot_grid(ggdraw() + draw_label(""),
                       ggdraw() + draw_label("Cell Ranger Metrics Summary"),
                       plot.tab(plotdf),
                       ggdraw() + draw_label(""),
                       ncol=1, rel_heights=c(0.1,0.1,0.7,0.2))
      make.pdf(plot, paste0(fig_path, "/0cellranger.pdf"), 7, 8)
    }
  }, error = function(e) {
    message("An error occurred in 0.cellranger output: ", e$message)
  })
  
  ###### Page 1: cell calling
  tryCatch({
    message('1.cell calling')
    if (file.exists(Molpath)) {
      plot = UvsI(obj, Molpath) 
    } else {
      plot = gdraw("No molecule_info.h5 found")
    }
    make.pdf(plot, paste0(fig_path, "/1cellcalling.pdf"), 7, 8)
    gc()
  }, error = function(e) {
    message("An error occurred in 1.cell calling: ", e$message)
  })
  
  ###### Page 2: UMAP + metrics 
  tryCatch({
    message('2.UMAP and metrics')
    plot = plot_grid(DimPlot(obj,label=T)+
                       ggtitle(g("UMAP"))+
                       NoLegend()+
                       theme(plot.title=element_text(hjust=0.5), 
                             axis.title.x=element_blank(), 
                             axis.title.y=element_blank())+
                       coord_fixed(ratio=1),
                     VlnPlot(obj,"logumi")+
                       NoLegend()+
                       theme(plot.title=element_text(hjust=0.5), 
                             axis.title.x=element_blank(), 
                             axis.title.y=element_blank()),
                     FeaturePlot(obj,"percent.mt")+
                       ggtitle("%MT")+
                       theme(plot.title=element_text(hjust=0.5), 
                             axis.title.x=element_blank(), 
                             axis.title.y=element_blank())+
                       coord_fixed(ratio=1)+
                       theme(legend.position="top",
                             legend.justification="center",
                             legend.key.width=unit(2, "lines")),
                     FeaturePlot(obj,"pct.intronic")+
                       ggtitle("%Intronic")+
                       theme(plot.title=element_text(hjust=0.5), 
                             axis.title.x=element_blank(), 
                             axis.title.y=element_blank())+
                       coord_fixed(ratio=1)+
                       theme(legend.position="top",
                             legend.justification="center",
                             legend.key.width=unit(2,"lines")),
                     ncol=2)
    make.pdf(plot, paste0(fig_path, "/2umap.pdf"), 7, 8)
  }, error = function(e) {
    message("An error occurred in 2.UMAP and metrics: ", e$message)
  })
  
  ###### Page 3: Raw spatial data
  tryCatch({
    message('3.Raw spatial data')
    plot <- make_spatial_barcode_plot(df, puckdf, SBh5path)
    make.pdf(plot, paste0(fig_path, "/3rawspatial.pdf"),7,8)
  }, error = function(e) {
    message("An error occurred in 3.Raw spatial data: ", e$message)
  })
  
  ###### Page 4: Beadplot
  tryCatch({
    message('4.Beadplot')
    sb.data = df %>% count_umis %>% 
      merge(y=puckdf, all.x=T, by="sb_index") %>% 
      group_by(sb_index) %>% 
      summarize(umi=sum(umi), x_um=mean(x_um), y_um=mean(y_um)) %>%
      ungroup %>% 
      filter(!is.na(x_um),!is.na(y_um)) %>% 
      arrange(umi)
    
    p1 = beadplot(sb.data, max(sb.data$umi), "raw")
    make.pdf(p1, paste0(fig_path, "/4beadplot.pdf"), 7, 8)
  }, error = function(e) {
    message("An error occurred in 4.Beadplot: ", e$message)
  })
  
  ###### Page 5: DBSCAN
  tryCatch({
    message('5.DBSCAN')
    plot <- plot_dbscan(obj, coords)
    make.pdf(plot, paste0(fig_path, "/5DBSCAN.pdf"),7,8)
  }, error = function(e) {
    message("An error occurred in 5.DBSCAN: ", e$message)
  })
  
  ###### Page 6: Spatial
  tryCatch({
    message('6.Spatial')
    p1 = DimPlot(obj[,!is.na(obj$x_um)], reduction="spatial")+
      coord_fixed(ratio=1)+
      ggtitle(g("%placed: {round(sum(!is.na(coords$x_um))/nrow(coords)*100,2)} ({sum(!is.na(obj$x_um))}/{ncol(obj)})")) + 
      NoLegend() + 
      xlab("x-position (\u00B5m)") + 
      ylab("y-position (\u00B5m)")
    p2 = DimPlot(obj[,!is.na(obj$x_um)], reduction="spatial",split.by="seurat_clusters",ncol=5) + 
      theme_void() + 
      coord_fixed(ratio=1) +
      NoLegend()
    plot = plot_grid(p1, p2, ncol=1, rel_heights=c(1,1))
    make.pdf(plot, paste0(fig_path, "/6spatial.pdf"),7,8)
  }, error = function(e) {
    message("An error occurred in 6.Spatial: ", e$message)
  })
  
  ###### Page 7: Create metrics plot
  tryCatch({
    message('7.Create metrics plot')
    res = metrics_plots(obj)
    plot <- res[[1]]
    Misc(obj,"SB_filtering") <- res[[2]]
    make.pdf(plot, paste0(fig_path, "/7metrics.pdf"),7,8)
  }, error = function(e) {
    message("An error occurred in 7.Metrics plot: ", e$message)
  })

  ###### Page 8: Sample bead plots
  tryCatch({
    message('8.Sample bead plots')
    plots <- sample_bead_plots(data.list, puckdf)
    make.pdf(plots, paste0(fig_path, "/8SampleBead.pdf"),7,7)
  }, error = function(e) {
    message("An error occurred in 8.Sample bead plots: ", e$message)
  })
  
  message("Done!")
  
}
