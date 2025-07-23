#' Sub-sampled CosMx non-small cell lung cancer (NSCLC) spatial transcriptomics dataset
#'
#' A Bruker CosMx spatial transcriptomics dataset containing a sub-region of a non-small cell lung cancer (NSCLC) sample. The data contains pathology annotations. The full dataset is available through the SubcellularSpatialData R/Bioconductor package.
#'
#' @format A SpatialExperiment object.
#' @docType data
#' @references Bhuva DD, Tan CW, Salim A, Marceaux C, Pickering MA, Chen J, Kharbanda M, Jin X, Liu N, Feher K, Putri G. Library size confounds biology in spatial transcriptomics data. Genome Biology. 2024 Apr 18;25(1):99.
#' 
#'  Bhuva DD, Tan CW, Marceaux C, Pickering M, Salim A, Chen J, Kharbanda M, Jin X, Liu N, Feher K, et al. Library size confounds biology in spatial transcriptomics data. 2024. Zenodo. https://doi.org/10.5281/zenodo.7959786.
#' 
#'  Bhuva DD: SubcellularSpatialData: annotated spatial transcriptomics datasets from 10x Xenium, NanoString CosMx and BGI STOmics. Bioconductor. 2024 https://doi.org/10.18129/B9.bioc.SubcellularSpatialData.
#' 
#'  Bhuva DD. Library size confounds biology in spatial transcriptomics. 2024. Zenodo. https://doi.org/10.5281/zenodo.10946961.
#' 
"lung5"

#' Pseudo-bulked single-cell non-small cell lung cancer (NSCLC) atlas reference dataset
#'
#' A pseudo-bulked single-cell RNA-seq non-small cell lung cancer (NSCLC) atlas dataset from the study by Salcher et al., Cancer Cell, 2022. The original dataset has been subsetted to the genes present in the 10x Xenium dataset.
#'
#' @format A SingleCellExperiment object.
#' @docType data
# 
#' @references Salcher S, Sturm G, Horvath L, Untergasser G, Kuempers C, Fotakis G, Panizzolo E, Martowicz A, Trebo M, Pall G, Gamerith G. High-resolution single-cell atlas reveals diversity and plasticity of tissue-resident neutrophils in non-small cell lung cancer. Cancer cell. 2022 Dec 12;40(12):1503-20.
#' 
"ref_luca"

.myDataEnv <- new.env(parent = emptyenv()) # not exported

.data_internal <- function(dataset) {
  if (!exists(dataset, envir = .myDataEnv)) {
    utils::data(list = c(dataset), envir = .myDataEnv)
  }
}
