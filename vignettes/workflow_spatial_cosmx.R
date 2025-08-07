## ----setup, include=FALSE-----------------------------------------------------
#set knitr chunk options
knitr::opts_chunk$set(warning = FALSE, message = FALSE)

#load packages to avoid startup messages later in the code
library(ggplot2)
library(patchwork)
library(scran)
library(scater)
library(scuttle)
library(standR)
library(SingleR)
library(SpatialExperiment)
library(SpaNorm)
library(singscore)
library(CosMxSpatialAnalysisWorkshop)

#automatically create a bib database for R packages
knitr::write_bib(c(
  .packages(), 'knitr', 'rmarkdown', 'prettydoc'
), 'packages.bib')

## -----------------------------------------------------------------------------
library(SpatialExperiment)

data(lung5) # replicate 1
# preview object
lung5
# retrieve column (cell) annotations
colData(lung5)
# retrieve row (gene/probe) annotations
rowData(lung5)
# retrieve spatial coordinates
head(spatialCoords(lung5))
# retrieve counts - top 5 cells and genes
counts(lung5)[1:5, 1:5]

## -----------------------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(SpaNorm)

# set default point size to 0.5 for this report
update_geom_defaults(geom = "point", new = list(size = 0.5))
plotSpatial(lung5, colour = region) +
  scale_colour_brewer(palette = "Dark2") +
  # improve legend visibility
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5)))

## ----eval = FALSE-------------------------------------------------------------
# ## DO NOT RUN during the workshop as the data are large and could crash your session
# if (!require("ExperimentHub", quietly = TRUE)) {
#   BiocManager::install("ExperimentHub")
# }
# if (!require("SubcellularSpatialData", quietly = TRUE)) {
#   BiocManager::install("SubcellularSpatialData")
# }
# 
# # load the metadata from the ExperimentHub
# eh = ExperimentHub()
# # download the transcript-level data
# tx = eh[["EH8232"]]
# # filter the Lung5_Rep1 sample
# tx = tx[tx$sample_id == "Lung5_Rep1", ]
# # summarise measurements for each cell
# lung5 = tx2spe(tx, bin = "cell")

## -----------------------------------------------------------------------------
table(rowData(lung5)$genetype)

## -----------------------------------------------------------------------------
library(scater)
library(scuttle)

# identify the negative control probes
is_neg = rowData(lung5)$genetype != "Gene"
# compute QC metrics for each cell
lung5 = addPerCellQCMetrics(lung5, subsets = list("neg" = is_neg))

## -----------------------------------------------------------------------------
# compute log counts - to fix the skew with counts
logcounts(lung5) = log2(counts(lung5) + 1)
# compute PCA
lung5 = runPCA(lung5, ncomponents = 30)
# compute UMAP
lung5 = runUMAP(lung5, dimred = "PCA")

## -----------------------------------------------------------------------------
plotColData(lung5, "detected", "sum", colour_by = "region", scattermore = TRUE) +
  scale_colour_brewer(palette = "Dark2") +
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5), title = "Region")) +
  geom_hline(yintercept = 50, lty = 2) +
  geom_vline(xintercept = 75, lty = 2)

## -----------------------------------------------------------------------------
library(standR)

# UMAP plot
p1 = plotDR(lung5, dimred = "UMAP", colour = sum, alpha = 0.75, size = 0.5) +
  scale_colour_viridis_c(option = "F") +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "bottom")
# spatial plot
p2 = plotSpatial(lung5, colour = sum, alpha = 0.75) +
  scale_colour_viridis_c(option = "F") +
  theme(legend.position = "bottom")
p1 + p2 + plot_annotation(title = "Sum (library size)")

## -----------------------------------------------------------------------------
# UMAP plot
p1 = plotDR(lung5, dimred = "UMAP", colour = detected, alpha = 0.75, size = 0.5) +
  scale_colour_viridis_c(option = "F") +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "bottom")
# spatial plot
p2 = plotSpatial(lung5, colour = detected, alpha = 0.75) +
  scale_colour_viridis_c(option = "F") +
  theme(legend.position = "bottom")
p1 + p2 + plot_annotation(title = "Detected (number of genes expressed)")

## -----------------------------------------------------------------------------
# UMAP plot
p1 = plotDR(lung5, dimred = "UMAP", colour = subsets_neg_percent, alpha = 0.75, size = 0.5) +
  scale_colour_viridis_c(option = "F", limits = c(0, 5), oob = scales::squish) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "bottom")
# spatial plot
p2 = plotSpatial(lung5, colour = subsets_neg_percent, alpha = 0.75) +
  scale_colour_viridis_c(option = "F", limits = c(0, 5), oob = scales::squish) +
  theme(legend.position = "bottom")
p1 + p2 + plot_annotation(title = "Proportion of negative control probes")

## -----------------------------------------------------------------------------
thresh_sum = 75
thresh_detected = 50
thresh_neg = 5

# classify cells to keep
lung5$keep = lung5$sum > thresh_sum & lung5$detected > thresh_detected & lung5$subsets_neg_percent < thresh_neg

# UMAP plot
p1 = plotDR(lung5, dimred = "UMAP", colour = keep, alpha = 0.75, size = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5)))
# spatial plot
p2 = plotSpatial(lung5, colour = keep, alpha = 0.75) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5)))
p1 + p2 + plot_annotation(title = "Cells retained after QC filtering")

## -----------------------------------------------------------------------------
# apply the filtering
lung5 = lung5[, lung5$keep]
# remove negative probes
lung5 = lung5[!is_neg, ]

## -----------------------------------------------------------------------------
lung5$logLS = log1p(lung5$sum) - mean(log1p(lung5$sum))
plotSpatial(lung5, colour = logLS) +
  scale_colour_distiller(palette = "PRGn", limits = c(-0.5, 0.5), oob = scales::squish) +
  theme(legend.position = "bottom")

