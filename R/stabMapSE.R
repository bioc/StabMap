#' Stabilised mosaic single cell data integration using unshared features
#'
#' stabMapSE performs StabMap with SummarizedExperiment type
#' (SingleCellExperiment, SpatialExperiment) objects as input.
#'
#' @param SE_list List of SummarizedExperiment, SingleCellExperiment,
#' SpatialExperiment, etc that should be named.
#' @param assays Named character vector of assays to be extracted from the
#' objects.
#'
#' @return matrix containing common embedding with rows corresponding to cells,
#' and columns corresponding to PCs or LDs for reference dataset(s).
#' @importFrom SummarizedExperiment assay
stabMapSE <- function(SE_list,
                      assays = "logcounts",
                      args = list()) {

  browser()

  if (length(assays) == 1) {
    assays <- setNames(rep(assays, length(SE_list)), names(SE_list))
  }

  stopifnot(
    "SE_assay_names should have names corrisponding to the SE objects in assay_list" = setequal(names(SE_list), names(assays))
    )

  assay_list <- mapply(
    SummarizedExperiment::assay,
    SE_list,
    assays[names(SE_list)],
    SIMPLIFY = FALSE
  )

  return(assay_list)
}
