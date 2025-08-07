## -----------------------------------------------------------------------------
# remotes::install_github("bhuvad/CosMxSpatialAnalysisWorkshop", ref = "PhiSpace")

# Key methods
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(SpaNorm))
# Plotting and publishing results
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(seriation))

suppressPackageStartupMessages(library(CosMxSpatialAnalysisWorkshop))

## Some parameters and functions
tissueNames_CosMx <- c(
  "Lung5_Rep1", "Lung5_Rep2", "Lung5_Rep3",
  "Lung6", "Lung9_Rep1", "Lung9_Rep2", "Lung12", "Lung13"
)

# Save PhiSpace heatmaps 
tempSavePlots <- function(
    spe, methodName = "PhiSpace", tissueName = "Lung", coordNames = c("x", "y"), 
    saveDir = "~/Desktop/spatialCellTypeHeatmaps/",
    legendPosition = "none", censQuant = 1, freeColScale = F, quants = c(0.1,1), psize = 0.5,
    width = 10, height = 10, fignrow = 4, figncol = 4, fsize_title = 10, plot_margin = c(0,0,0,0)
){
  
  scores <- as.matrix(reducedDim(spe, methodName))
  
  if(freeColScale){
    lmts <- NULL
  } else {
    lmts <- quantile(scores, quants)
  }
  
  ctypes <- colnames(scores) %>% sort()
  outPlots <- vector("list", length(ctypes)) 
  names(outPlots) <- ctypes
  for(i in 1:length(ctypes)){
    ctype <- ctypes[i]
    cols <- scores[, ctype]
    cols <- PhiSpace:::censor(cols, quant = censQuant)
    df <- spe@colData %>% as.data.frame()
    df <- cbind(df, spatialCoords(spe))
    outPlots[[i]] <- df %>% mutate(cols = cols) %>% arrange(cols) %>%
      ggplot(aes(x = !!sym(coordNames[1]), y = !!sym(coordNames[2]))) +
      geom_point(aes(colour = cols), size = psize, alpha = 1, stroke = 0) +
      theme_void() + scale_colour_gradientn(colours = PhiSpace:::MATLAB_cols, limits = lmts) + 
      ggtitle(ctype) + theme(
        legend.position = legendPosition, 
        plot.title = element_text(size = fsize_title, hjust = 0.5, vjust = -1),
        plot.margin = unit(plot_margin, "points")
      )
  }
  
  npheno <- length(outPlots)
  nFigPerPlot <- fignrow * figncol
  for(ii in 1:ceiling(npheno/nFigPerPlot)){
    idx_start <- (ii-1)*nFigPerPlot + 1
    idx_end <- min(ii*nFigPerPlot, npheno)
    ggsave(
      paste0(saveDir, tissueName, "/", tissueName, "_", methodName, "_", ii, ".png"),
      ggpubr::ggarrange(plotlist = outPlots[idx_start:idx_end], ncol = figncol, nrow = fignrow), 
      width = width, height = height
    )
  }
}

# Plot cell type copresence matrix as heatmap
tempCopresence <- function(interaction_matrix, fsize = 12, useSeriate = T){
  o <- seriate(interaction_matrix)
  hm <- Heatmap(
    interaction_matrix, name = "Correlation",
    show_row_names = T, show_column_names = F, show_heatmap_legend = F,
    row_order = get_order(o, 1), column_order = get_order(o, 2), 
    row_names_gp = gpar(fontsize = fsize), column_names_gp = gpar(fontsize = fsize)
  )
  
  return(hm)
}

## -----------------------------------------------------------------------------
data("lung5_norm")
data("ref_luca")
# load("~/Desktop/MIG_spatial_data/ref_luca.rda")
# load("~/Desktop/MIG_spatial_data/lung5_norm.rda")

## -----------------------------------------------------------------------------
ref_luca <- scranTransf(ref_luca)
# lung5_norm already normalised by SpaNorm

## -----------------------------------------------------------------------------
colData(ref_luca) 
YtrainName <- "cell_type"

## -----------------------------------------------------------------------------
lung5_norm <- PhiSpace(
  reference = ref_luca, query = lung5_norm,
  phenotypes = YtrainName, refAssay = "logcounts", regMethod = "PLS"
)
reducedDim(lung5_norm, "PhiSpace")[1:5,1:5]

## -----------------------------------------------------------------------------
lung5_norm$PhiCellType <- getClass(reducedDim(lung5_norm, "PhiSpace"))

## -----------------------------------------------------------------------------
p <- VizSpatial(lung5_norm, ptSize = 1, groupBy = "PhiCellType") 
p
plotly::ggplotly(p)

## ----eval=FALSE---------------------------------------------------------------
# tempSavePlots(lung5_norm, tissueName = "lung5_norm_sub", freeColScale = T, censQuant = 0.5, psize = 0.8)

## ----eval=FALSE---------------------------------------------------------------
# # sc_lung5_norm_sub <- readRDS("~/Desktop/MIG_spatial_data/lung5_norm_sub_PhiSpace.rds")
# # use_data(sc_lung5_norm_sub, overwrite = T)
# # load("~/Desktop/MIG_spatial_data/sc_lung5_norm_sub.rda")
# data("sc_lung5_norm_sub")
# reducedDim(lung5_norm, "Phi5Ref") <- sc_lung5_norm_sub
# tempSavePlots(lung5_norm, tissueName = "lung5_norm_sub", methodName = "Phi5Ref", freeColScale = T, censQuant = 0.5, psize = 0.8)

