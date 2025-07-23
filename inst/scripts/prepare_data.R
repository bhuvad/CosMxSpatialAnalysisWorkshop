library(ggplot2)
library(patchwork)
library(scran)
library(scater)
library(scuttle)
library(standR)
library(SpatialExperiment)
library(SpaNorm)
library(singscore)
library(ExperimentHub)
library(SubcellularSpatialData)
library(BiocParallel)

# load transcript data
eh = ExperimentHub()
tx = eh[["EH8232"]]
tx$fov = paste0("fov", tx$fov)
spe = tx2spe(tx, bin = "cell")
saveRDS(spe, file = "cosmx_hs_nsclc_spe.rds")

# process all data
speapply <- \(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, foreach.args = list()) {
  require("foreach")
  stopifnot(is(X, "SingleCellExperiment"))

  # split
  sample_ids <- unique(X$sample_id)
  if (USE.NAMES) {
    names(sample_ids) <- sample_ids
  }
  X <- lapply(sample_ids, \(sid) X[, X$sample_id == sid])

  # apply
  foreach.args <- foreach.args[setdiff(names(foreach.args), c("x", ".combine"))]
  foreach.args <- c(list(x = iterators::iter(X)), foreach.args)
  X <- do.call(what = foreach::foreach, args = foreach.args) %dopar% {
    FUN(x, ...)
  }
  names(X) <- names(sample_ids)

  # merge
  if (simplify & all(sapply(X, is, "SingleCellExperiment"))) {
    X <- do.call(cbind, X)
  } else if (simplify | simplify == "array") {
    higher <- simplify == "array"
    X <- simplify2array(X, higher = higher)
  }

  X
}

# QC
thresh_sum = 50
thresh_detected = 20
thresh_neg = 5
is_neg = rowData(spe)$genetype != "Gene"
spe = addPerCellQCMetrics(spe, subsets = list("neg" = is_neg))
spe$keep = spe$sum > thresh_sum & spe$detected > thresh_detected & spe$subsets_neg_percent < thresh_neg
spe = spe[!is_neg, spe$keep]

# normalisation
speapply(spe, \(x) {
  cl = scran::quickCluster(x, BPPARAM = MulticoreParam(4))
  x = computeSumFactors(x, clusters = cl, BPPARAM = MulticoreParam(4))
  set.seed(1000)
  # perform SpaNorm normalisation
  x = SpaNorm(x, df.tps = c(10, 3), sample.p = 0.05)
  saveRDS(x, sprintf("~/Downloads/%s.rds", x$sample_id[[1]]))
  gc()
})
