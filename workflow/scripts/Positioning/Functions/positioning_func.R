#########################################################################################
################################ POSITIONING FUNCTIONS ##################################
#########################################################################################

# helper methods
chunk_vector <- function(v, chunk_size) {return(split(v, ceiling(seq_along(v) / chunk_size)))}
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
h_index <- function(vec) {
  hdf = data.frame(x=vec) %>% group_by(x) %>% summarize(n=n()) %>% arrange(desc(x)) %>% mutate(n = cumsum(n))
  return(max(pmin(hdf$x,hdf$n)))
}


# Do a grid search to find the ideal DBSCAN parameters
opt_dbscan <- function(data.list) {
  res = data.frame()
  for (k in 0:26) {
    params = expand.grid(eps.vec, minPts.vec) %>% setNames(c("eps","minPts"))
    row_lists = chunk_vector(1:nrow(params), ceiling(nrow(params)/ncores))
    
    params$pct = furrr::future_map(row_lists, function(v) {
      map_dbl(v, function(i) {
        m = map_int(data.list, ~dbscan::dbscan(.[c("x_um","y_um")], eps=params$eps[[i]], minPts=params$minPts[[i]], weights=.$umi)$cluster %>% max)
        return(sum(m==1)/length(m))
      })
    }, .options=furrr_options(seed=T)) %>% flatten_dbl
    
    res = rbind(res, params)
    
    p = res %>% group_by(minPts) %>% summarize(pct=max(pct), .groups="drop") %>% arrange(minPts) %>% pull(pct) %>% {which.max(.)/len(.)}
    if (p < 0.9) {
      break
    }
    minPts.vec %<>% add(len(minPts.vec))
  }
  res$is.max = res$pct==max(res$pct)
  return(res)
}



# Assign centroid and record metadata
create_dbscan_coords <- function(data.list) {
  coords = lapply(data.list, function(df) {
    p = c(x_um = NA,
          y_um = NA,
          clusters = max(df$cluster),
          umi = sum(df$umi),
          beads = nrow(df),
          SNR = NA,
          minPts = unique(df$minPts) %>% {ifelse(is.null(.), NA, .)} %T>% {stopifnot(len(.)==1)},
          eps = unique(df$eps) %>% {ifelse(is.null(.), NA, .)} %T>% {stopifnot(len(.)==1)},
          pct.placed = unique(df$pct.placed) %>% {ifelse(is.null(.), NA, .)} %T>% {stopifnot(len(.)==1)}
         )
    if (max(df$cluster) == 1) {
      sdf = dplyr::filter(df, cluster==1)
      p[["x_um"]] = weighted.mean(sdf$x_um, w=sdf$umi)
      p[["y_um"]] = weighted.mean(sdf$y_um, w=sdf$umi)
      p[["SNR"]] = sum(sdf$umi)/sum(df$umi)
    }
    return(p)
  }) %>% bind_rows %>% as.data.frame %>% mutate(cb_index=as.numeric(names(data.list))) %>% select(cb_index, everything())
  return(coords)
}



# Plot the results
plot_dbscan <- function(coords, optim_plot) {
  # Panel 1: DBSCAN cluster distribution
  d = data.frame(x=coords$clusters) %>% rowwise %>% mutate(x=min(x,5)) %>% ungroup
  p1 = ggplot(d, aes(x=x)) + geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))*100), binwidth=.5) +
    geom_text(aes(label = sprintf("%1.0f%%", after_stat(count)/sum(after_stat(count))*100), y=after_stat(count)/sum(after_stat(count))*100), stat="bin", binwidth=1, vjust=-0.5)+
    theme_classic() + xlab("DBSCAN clusters") + ylab("Percent") +
    scale_y_continuous(limits=c(0,100)) +
    scale_x_continuous(breaks=min(d$x):max(d$x), labels=(min(d$x):max(d$x)) %>% {ifelse(.==5, "5+", .)}) +
    ggtitle("Cluster distribution")
  
  # Panel 3: SB UMI distribution
  d = coords %>% rowwise %>% mutate(x=min(clusters,5)) %>% ungroup
  p3 = ggplot(d, aes(x=as.factor(x), y=log10(umi))) + geom_violin(scale="count") + 
    scale_x_discrete(breaks=min(d$x):max(d$x), labels=(min(d$x):max(d$x)) %>% {ifelse(.==5, "5+", .)}) +
    xlab("DBSCAN clusters") + ylab("log10 SB UMI") + ggtitle("SB UMI distribution") + theme_classic()
  
  # Panel 4: SNR density
  max_density_x = mean(coords$SNR, na.rm=T)
  p4 = coords %>% filter(!is.na(x_um)) %>% ggplot(aes(x = SNR)) + geom_density() + 
    theme_minimal() +
    labs(title = "SNR per cell (density)", x = "SNR", y = "Density") + 
    geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
    annotate(geom = 'text', label = g("Mean: {round(max_density_x*100, 1)}%"), x = max_density_x+0.01, y = Inf, hjust = 0, vjust = 1, col="red")

  plot = plot_grid(gdraw("DBSCAN Results"),
                   plot_grid(p1, optim_plot, p3, p4, ncol=2),
                   ncol=1, rel_heights=c(0.05,0.95))
  return(plot)
}



# note: density of lone point is num umi, then add umis from surroundings weighted by umi
# note: filter out singletons when max umi > 2
kde <- function(df, bw, radius) {
  cb_i = unique(df$cb_index) ; stopifnot(len(cb_i)==1)
  df %<>% select(x_um, y_um, umi)
  res = c(cb_index=cb_i,
          x_um=NA, y_um=NA, x2_um=NA, y2_um=NA,
          d1=NA, d2=NA,
          sumi=NA, sbeads=NA, max=max(df$umi), r=NA, h=NA)
  if(nrow(df) == 1) {
    res[c("x_um", "y_um", "d1", "d2", "umi", "beads", "r", "h")] = c(df$x_um, df$y_um, df$umi, 0, df$umi, 1, 1, 1)
    return(res)
  }
  if (max(df$umi) >= 3) {df %<>% filter(umi > 1)}
  xmu = pdist(df[,c("x_um","y_um")])
  df$density = exp(-xmu^2/bw) %>% sweep(MARGIN=1, STATS=df$umi, FUN="*") %>% colSums # broadcast the umi vector across the columns
  
  rowmax = df %>% {.[which.max(.$density),]}
  df$near = (df$x_um-rowmax$x_um)^2 + (df$y_um-rowmax$y_um)^2 < radius^2
  rowmax2 = df %>% filter(!near) %>% {.[which.max(.$density),]}
  if (nrow(rowmax2)==0) {df$near2=F} else {df$near2 = (df$x_um-rowmax2$x_um)^2 + (df$y_um-rowmax2$y_um)^2 < radius^2}
  
  # if (any(df$near & df$near2)) {print("check debug, might increase radius")}
  
  sdf <- df %>% filter(near)
  sdf2 <- df %>% filter(near2)
  
  res[c("x_um","y_um","x2_um","y2_um")] = c(weighted.mean(sdf$x_um, w=sdf$umi),
                                            weighted.mean(sdf$y_um, w=sdf$umi),
                                            weighted.mean(sdf2$x_um, w=sdf2$umi) %>% {ifelse(is.nan(.), NA, .)},
                                            weighted.mean(sdf2$y_um, w=sdf2$umi) %>% {ifelse(is.nan(.), NA, .)}
  )
  res[c("d1","d2","sumi","sbeads","r","h")] = c(rowmax$density,
                                                rowmax2$density %>% {ifelse(len(.)==0, 0, max(.))},
                                                sum(sdf$umi),
                                                nrow(sdf),
                                                r=sum(sdf$umi)/max(sdf$umi),
                                                h=h_index(sdf$umi)
  )
  return(res)
}



plot_kde <- function(kde_coords) {
  p1 <- kde_coords %>% ggplot(aes(x=x_um,y=y_um))+geom_point(size=0.1,shape=16)+coord_fixed()+theme_classic()+ggtitle("Location of highest density")+xlab("")+ylab("")+theme(plot.title=element_text(size=12))
  p3 <- kde_coords %>% filter(!is.na(x2_um),!is.na(y2_um)) %>% ggplot(aes(x=x2_um,y=y2_um))+geom_point(size=0.1,shape=16)+coord_fixed()+theme_classic()+ggtitle("Location of second-highest density")+xlab("")+ylab("")+theme(plot.title=element_text(size=12))
  #p1 <- kde_coords %>% ggplot(aes(x=x_um,y=y_um))+geom_bin2d(bins=c(round(xrange/50),round(yrange/50)))+coord_fixed()+theme_classic()+ggtitle("Distribution of highest density")+xlab("")+ylab("")+theme(plot.title=element_text(size=12))
  #p3 <- kde_coords %>% ggplot(aes(x=x2_um,y=y2_um))+geom_bin2d(bins=c(round(xrange/50),round(yrange/50)))+coord_fixed()+theme_classic()+ggtitle("Distribution of second-highest density")+xlab("")+ylab("")+theme(plot.title=element_text(size=12))
  
  label = g("3:1\ncells: {sum(kde_coords$ratio<1/3)}/{nrow(kde_coords)} ({round(sum(kde_coords$ratio<1/3)/nrow(kde_coords)*100,2) %>% paste0('%')})")
  p2 <- kde_coords %>% ggplot(aes(x = ratio)) + geom_density() + 
    theme_minimal() + theme(plot.title=element_text(size=12)) + 
    labs(title = "Distribution of top-2 density ratio", x = "Ratio", y = "Density") + 
    geom_vline(xintercept = 1/3, color = "red", linetype = "dashed") +
    annotate(geom = 'text', label = label, x = 1/3+0.01, y = Inf, hjust = 0, vjust = 1, col="red")
  
  d <- kde_coords %>% filter(!is.na(x2_um),!is.na(y2_um)) %>% mutate(dist=sqrt((x_um-x2_um)^2+(y_um-y2_um)^2))
  p4 <- ggplot(d, aes(x=dist)) + geom_histogram(bins = 30) + theme_bw() +
    xlab("Distance (\u00B5m)")+ylab("Frequency")+ggtitle("Distance between top-2 densities")+theme(plot.title=element_text(size=12))
  
  plot = plot_grid(gdraw("KDE Results"),
                   plot_grid(p1, p2, p3, p4, ncol=2),
                   ncol=1, rel_heights=c(0.05,0.95))
  return(plot)
}



dbscan_vs_kde <- function(coords) {
  coords %<>% mutate(dist=sqrt((x_um_dbscan-x_um_kde)^2+(y_um_dbscan-y_um_kde)^2))
  
  # Panel 1: Distance between DBSCAN and KDE assignments
  d = filter(coords, !is.na(dist))
  max_density_x = median(d$dist+1, na.rm=T) %>% log10
  p1 <- ggplot(d, aes(x=log10(dist+1)))+geom_histogram(bins=30)+theme_bw()+
    labs(title = "DBSCAN vs KDE distance", x = "log1p Distance (\u00B5m)", y = "Frequency") + 
    geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
    annotate(geom = 'text', label = round(10^max_density_x-1, 2) %>% paste0("\u00B5m"), x = max_density_x+0.1, y = Inf, hjust = 0, vjust = 1.3, col="red")
  
  # Panel 2: Distribution of KDE ratio for each DBSCAN cluster
  d = coords %>% rowwise %>% mutate(clusters=min(clusters,5)) %>% ungroup
  p2 <- ggplot(d, aes(x=clusters %>% as.factor, y=ratio)) + geom_violin() +
    theme_classic() + xlab("DBSCAN clusters") + ylab("KDE ratio") +
    scale_x_discrete(breaks=min(d$clusters):max(d$clusters), labels=(min(d$clusters):max(d$clusters)) %>% {ifelse(.==5, "5+", .)}) +
    ggtitle("KDE ratio per DBSCAN cluster") + 
    geom_hline(yintercept = 1/3, color = "red", linetype = "dashed")
  
  # Panel 3: KDE ratio for disagreeing placements
  p3 <- coords %>% filter(clusters==1) %>% ggplot(aes(x=dist, y=ratio))+geom_point(size=0.5)+theme_bw()+
    xlab("Distance between assignments")+ylab("KDE ratio")+ggtitle("Density ratio of disagreements")+
    geom_hline(yintercept = 1/3, color = "red", linetype = "dashed")
  
  # Panel 4: Contingency table of placements
  d = coords %>% mutate(dbscan_pass=clusters==1, kde_pass=ratio<1/3) %>% group_by(dbscan_pass, kde_pass) %>% summarize(n=n(), .groups="drop") %>% mutate(pct=g("{round(n/sum(n)*100,2)}%\n{n}"))
  p4 <- ggplot(d, aes(x=dbscan_pass,y=kde_pass,fill=n))+geom_tile()+geom_text(label=d$pct)+theme_bw()+
    xlab("DBSCAN=1")+ylab("KDE < 1/3")+theme(legend.position="none")+coord_fixed()+ggtitle("Placement table")
    
  plot = plot_grid(gdraw("DBSCAN vs. KDE Comparison"),
                   plot_grid(p1, p2, p3, p4, ncol=2),
                   ncol=1, rel_heights=c(0.05,0.95))
  return(plot)
}


sample_bead_plots <- function(data.list, coords) {
  plot.sb <- function(cb_i) {
    row = coords %>% filter(cb_index==cb_i) %T>% {stopifnot(nrow(.)==1)}
    subdf = data.list[[as.character(cb_i)]] ; stopifnot(subdf$cb_index == cb_i)
    subdf %<>% arrange(umi) %>% filter(!is.na(x_um), !is.na(y_um))
    subdf1 <- filter(subdf, cluster==1)
    subdf2 <- filter(subdf, cluster==2)
    
    plot <- ggplot() + coord_fixed(ratio=1 ,xlim=xlims, ylim=ylims) + theme_void() +
      geom_point(data=subdf, mapping=aes(x=x_um, y=y_um, col=umi), size=2, shape=16)+
      geom_point(aes(x=row$x_um_kde, y=row$y_um_kde), color="red", shape=5, size=3) + 
      ggtitle(g("[{unique(subdf$cb_index)}] ({round(row$ratio,2)})")) +
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
            legend.box.spacing=unit(0,"cm"),
            plot.title = element_text(hjust = 0.5))
    # Add the DBSCAN=1 point
    if(nrow(subdf1) > 0) {
      plot <- plot + geom_point(aes(x=weighted.mean(subdf1$x_um, w=subdf1$umi),
                                    y=weighted.mean(subdf1$y_um, w=subdf1$umi)),
                                    color="red", shape=0, size=3)
    }
    # Add the DBSCAN=2 point
    if(nrow(subdf2) > 0) {
      plot <- plot + geom_point(aes(x=weighted.mean(subdf2$x_um, w=subdf2$umi),
                                    y=weighted.mean(subdf2$y_um, w=subdf2$umi)),
                                    color="green", shape=0, size=3)
    }
    return(plot)
  }
  
  list0 = data.list %>% keep(~max(.$cluster)==0 & nrow(.)>0) %>% names %>% as.numeric
  list1 = data.list %>% keep(~max(.$cluster)==1 & nrow(.)>0) %>% names %>% as.numeric
  list2 = data.list %>% keep(~max(.$cluster)==2 & nrow(.)>0) %>% names %>% as.numeric

  if(len(list0) > 0) {
    p1 = plot_grid(ggdraw()+draw_label("DBSCAN=0"),
                   list0 %>% sample(min(12,len(list0)), replace=F) %>% map(plot.sb) %>% {plot_grid(plotlist=.,ncol=3)}, ncol=1, rel_heights=c(0.1,2))
  } else {p1 = gdraw("No DBSCAN = 0")}
  if(len(list1) > 0) {
    p2 = plot_grid(ggdraw()+draw_label("DBSCAN=1"),
                   list1 %>% sample(min(12,len(list1)), replace=F) %>% map(plot.sb) %>% {plot_grid(plotlist=.,ncol=3)}, ncol=1, rel_heights=c(0.1,2))
  } else {p2 = gdraw("No DBSCAN = 1")}
  if(len(list2) > 0) {
    p3 = plot_grid(ggdraw()+draw_label("DBSCAN=2"),
                   list2 %>% sample(min(12,len(list2)), replace=F) %>% map(plot.sb) %>% {plot_grid(plotlist=.,ncol=3)}, ncol=1, rel_heights=c(0.1,2))
  } else {p3 = gdraw("No DBSCAN = 2")}
  
  plots = list(p1, p2, p3)
  
  return(plots)
}
