# set Lib Path
lib_path <- Sys.getenv('R_LIBS', unset = NA)
if (is.na(lib_path)) {stop("R_LIBS environment variable is not set.")}
.libPaths(lib_path)

# load packages
# Run positioning for each cell in input dataframe (matrix.csv -> coords.csv)
suppressMessages(library(glue)) ; g=glue ; len=length
suppressMessages(library(gridExtra))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dbscan))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(Matrix))
suppressMessages(library(rdist))
suppressMessages(library(furrr))
suppressMessages(library(future))
suppressMessages(library(parallel))
options(warn = -1)


# Get the path to the R function
r_func_path <- Sys.getenv('R_FUNC', unset = NA)
if (is.na(r_func_path)) {stop("R_FUNC environment variable is not set.")}
source(paste0(r_func_path, '/Functions/positioning_func.R'))


# DBSCAN hyperparameters
eps.vec = c(50, 100)
minPts.vec = c(3:42) # auto-searches up to 40x26
# KDE hyperparameters
bw = 800     # bandwith of kernel
radius = 200 # inclusion radius around density peak

# parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 3) {
  matrix_path <- args[[1]] # (cell, x, y, umi) dataframe
  out_path <- args[[2]]    # path to write coords and plots
  ncores <- as.numeric(args[[3]])
} else {
  stop("Usage: Rscript positioning.R matrix_path output_path", call. = FALSE)
}
stopifnot(file.exists(matrix_path))
if (!dir.exists(out_path)) { dir.create(out_path, recursive = T) }


# load data.list
df = read.table(matrix_path, header=T, sep=",")
stopifnot(names(df) == c("cb_index", "x_um", "y_um", "umi"))
xlims = range(df$x_um) ; xrange = max(df$x_um) - min(df$x_um)
ylims = range(df$y_um) ; yrange = max(df$y_um) - min(df$y_um)
num_pucks = round(xrange/yrange) # for plot layout
data.list = split(df, df$cb_index)
data.list %<>% map(~arrange(.,desc(umi),cb_index))
rm(df) ; invisible(gc())
print(g("Running positioning on {len(data.list)} cells"))

# prepare for positioning
plan(multisession, workers=ncores)
params = opt_dbscan(data.list)
eps = params$eps[params$is.max][[1]]
minPts = params$minPts[params$is.max][[1]]
pct.placed = round(max(params$pct)*100, 2)
print(g("Optimal eps: {eps} \t Optimal minPts: {minPts} \t %placed: {pct.placed}"))
optim_plot = ggplot(params, aes(x=minPts, y=pct*100, col=as.factor(eps))) + geom_line() +
  theme_bw() + ylab("% Placed") + labs(col="eps") + ggtitle("Parameter optimization") +
  geom_vline(xintercept = minPts, color = "red", linetype = "dashed") +
  annotate(geom = 'text', label = g("eps: {eps}\nminPts: {minPts}\nplaced: {pct.placed}%"), x = minPts+1, y = 0, hjust = 0, vjust = 0, col="red") + 
  theme(legend.position = c(0.95, 0.50), legend.justification = c("right", "center"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"))
data.list %<>% lapply(function(df){
  mutate(df,
         cluster = dbscan::dbscan(df[c("x_um","y_um")], eps=eps, minPts=minPts, weights=df$umi, borderPoints=F)$cluster,
         eps = eps,
         minPts = minPts,
         pct.placed = pct.placed)
})

dbscan_coords <- create_dbscan_coords(data.list)
invisible(gc())

plot <- plot_dbscan(dbscan_coords, optim_plot)
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path, "DBSCAN.pdf"), 7, 8)))


### Run KDE ####################################################################
kde_coords <- map(data.list, ~kde(., bw, radius)) %>% bind_rows
stopifnot(!is.na(kde_coords$d1), !is.na(kde_coords$d2))
kde_coords %<>% mutate(ratio = d2/d1) %>% select(1:7, ratio, everything())
plot <- plot_kde(kde_coords)
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path, "KDE.pdf"), 7, 8)))


### More plots + save output ###################################################
stopifnot(dbscan_coords$cb_index == kde_coords$cb_index)
coords <- merge(dbscan_coords, kde_coords, by="cb_index", suffix=c("_dbscan","_kde"))
coords %<>% rename(x2_um_kde=x2_um, y2_um_kde=y2_um)
plot <- dbscan_vs_kde(coords)
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path, "DBSCANvsKDE.pdf"), 7, 8)))

plots <- sample_bead_plots(data.list, coords)
suppressMessages(suppressWarnings(make.pdf(plots, file.path(out_path, "beadplots.pdf"), 7, 8)))


# set x_um, y_um to be the filtered DBSCAN placements where ratio<1/3
coords %<>% mutate(x_um = ifelse(ratio<1/3, x_um_dbscan, NA),
                  y_um = ifelse(ratio<1/3, y_um_dbscan, NA)) %>% 
                  select(cb_index, x_um, y_um, everything())
write.table(coords, file.path(out_path, "coords.csv"), sep=",", row.names=F, col.names=T, quote=F)

### Position debugging (optional) ##############################################

# plot.sb <- function(subdf) {
#   subdf %<>% arrange(umi)
#   ggplot() + coord_fixed(ratio=1, xlim=xlims, ylim=ylims) + theme_void() +
#     geom_point(data=subdf, mapping=aes(x=x_um, y=y_um, col=umi), size=2, shape=20)
# }
# plot.kde <- function(subdf) {
#   if(nrow(subdf)==0) {return(gdraw("No points"))}
#   p = Nebulosa:::wkde2d(x=subdf$x_um, y=subdf$y_um, w=subdf$umi, n=200, lims=c(xlims, ylims)) %>% {transmute(reshape2::melt(as.matrix(.[[3]])), x_um=.[[1]][Var1], y_um=.[[2]][Var2], value=value)}
#   rowmax = p[which.max(p$value),]
#   ggplot(p, aes(x=x_um, y=y_um, fill=value))+geom_tile()+coord_fixed(ratio=1) +
#     annotate("path", x=rowmax$x_um+radius*cos(seq(0,2*pi,length.out=100)), y=rowmax$y_um+radius*sin(seq(0,2*pi,length.out=100)))
# }
# plot.metadata <- function(row) {
#   row %<>% select(-x_um, -y_um)
#   row1 <- select(row, 2:10)
#   d1 <- data.frame(names(row1), round(as.numeric(row1[1,]),2)) %>% setNames(c("Data","Value"))
#   row2 <- select(row, 11:22)
#   d2 <- data.frame(names(row2), round(as.numeric(row2[1,]),2)) %>% setNames(c("Data","Value"))
#   
#   plot_grid(gdraw(g("[{row$cb_index}]")),
#             plot_grid(plot.tab(d1), plot.tab(d2), ncol=2),
#             ncol=1, rel_heights=c(0.05,0.95))
# }
# debug_coords <- function(data.list, coords) {
#   library(shiny)
#   ui <- fluidPage(
#     fluidRow(
#       column(6, plotOutput("plot1", click = "plot1_click")),
#       column(6, plotOutput("plot2"))
#     ),
#     fluidRow(
#       column(6, plotOutput("plot3")),
#       column(6, plotOutput("plot4"))
#     )
#   )
#   
#   server <- function(input, output) {
#     output$plot1 <- renderPlot({ ggplot(coords, aes(x=x_um, y=y_um)) + geom_point() + coord_fixed() + theme_void() })
#     output$plot2 <- renderPlot({ plot.new() })
#     output$plot3 <- renderPlot({ plot.new() })
#     output$plot4 <- renderPlot({ plot.new() })
#     
#     observeEvent(input$plot1_click, {
#       click <- input$plot1_click
#       if(!is.null(click)) {
#         dists <- sqrt((coords$x_um - click$x)^2 + (coords$y_um - click$y)^2)
#         row <- coords[which.min(dists),]
#         cb_index <- row$cb_index %>% as.character
#         df <- data.list[[cb_index]]
#         output$plot2 <- renderPlot({ plot.sb(df) })
#         output$plot3 <- renderPlot({ plot.kde(df) })
#         output$plot4 <- renderPlot({ plot.metadata(row) })
#       }
#     })
#   }
#   print(shinyApp(ui = ui, server = server))
#   return(T)
# }
# debug_coords(data.list, coords)
