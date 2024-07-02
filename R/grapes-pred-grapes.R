#' Binary operator for model predictions on data
#'
#' This function performs model predictions via the \code{predict} function
#' for each column of data.
#'
#' @usage data \%pred\% models
#' @param data is a matrix with rows corresponding to features, and columns
#' corresponding to cells/observations
#' @param models is a list of univariate outcome models with the features as
#' explanatory variables
#'
#' @return a matrix with rows equal to \code{length(models)} and columns
#' corresponding to cells/observations
#'
#' @keywords internal
"%pred%" <- function(data, models) {
  do.call(rbind, lapply(models, stats::predict, t(data)))
}
