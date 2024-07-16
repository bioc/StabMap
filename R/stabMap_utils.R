#' selectFeatures
#'
#' For a given assay and set of features, perform variance ranking and select
#' a subset of features
#'
#' @param assay An assay matrix rows are features, columns are cells
#' @param features Character vector of the current features that are selected
#' @param maxFeatures Integer of the number of maxFeatures to select
#'
#' @return A character vector of the selected features according to variance
#' ranking.
#'
#' @keywords internal
selectFeatures = function(assay, features, maxFeatures) {

  if (!requireNamespace("scran", quietly = TRUE)) {
    stop("Install 'scran' select features when maxFeatures is selected")
  }

  genevars <- scran::modelGeneVar(
    assay[features, ],
    min.mean = 0
  )
  genevars_sorted <- genevars[
    order(genevars$bio, decreasing = TRUE),
  ]
  features_selected <- rownames(
    genevars_sorted
  )[seq_len(maxFeatures)]
  return(features_selected)
}
