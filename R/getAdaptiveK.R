#' Adaptive k selection for KNN classification
#'
#' Given an error matrix, identify the k that maximises the accuracy for cells
#' belonging to a provided labelling/grouping. If no labelling given, expect a
#' cell-cell similarity network to identify the k that maximises the accuracy
#' for cells within that neighbourhood. If neither are given, simply treat all
#' cells as if they have the same labelling/grouping
#'
#' @param E An error matrix with rows corresponding to cells and columns
#' corresponding to candidate k values, with values themselves corresponding to
#' error values (either binary for single classification, or continuous after
#' multiple classification).
#' @param labels Group labels for cells.
#' @param local A neighbourhood index representation, as typically output using
#' BiocNeighbors::findKNN().
#' @param outputPerCell Logical whether to return adaptive k for each cell, not
#' just for each label type (used for when labels is given).
#' @param ... Includes return_colnames, whether to give the colnames of the best
#' selected, or just the index, which is default TRUE.
#'
#' @return Vector of adaptive k values.
#'
#' @examples
#' E = matrix(runif(100),20,5)
#' colnames(E) <- paste0("K_", 1:5)
#'
#' # generate cell labels
#' labels = factor(rep(letters[1:2], each = 10))
#'
#' # generate nearest neighbourhood index representation
#' data = matrix(rpois(10*20, 10), 10, 20) # 10 genes, 20 cells
#' local = BiocNeighbors::findKNN(t(data), k = 5, get.distance = FALSE)$index
#'
#' best_k_labels = getAdaptiveK(E,
#'                              labels = labels)
#' best_k_local = getAdaptiveK(E,
#'                             local = local
#' )
#'
#' @export
getAdaptiveK = function(E,
                        labels = NULL,
                        local = NULL,
                        outputPerCell = TRUE,
                        ...) {

  # adaptive k selection for KNN classification
  # Given an error matrix E, with rows corresponding to cells
  # and columns corresponding to candidate k values, with values
  # themselves corresponding to error values (either binary
  # for single classification, or continuous after multiple
  # classification)
  # and given an optional factor labelling/grouping of cells
  # identify the k that maximises the accuracy for cells belonging
  # to that label/group
  # if no labelling given, expect a cell-cell similarity network
  # to identify the k that maximises the accuracy for cells within
  # that neighbourhood
  # if neither are given, simply treat all cells as if they have
  # the same labelling/grouping.

  # ... includes return_colnames, whether to give the
  # colnames of the best selected, or just the index,
  # which is default TRUE

  # if outputPerCell then return a vector of adaptive k
  # values for each cell, not just for each label type
  # (used for when labels is given)

  # if both labels and local given, labels will be
  # prioritised

  # local is a neighbourhood index representation
  # as typically output using BiocNeighbors::findKNN()

  # example data generation
  # data = matrix(rpois(10*20, 10), 10, 20) # 10 genes, 20 cells
  # local = BiocNeighbors::findKNN(t(data), k = 5, get.distance = FALSE)$index
  # E = matrix(runif(100),20,5)
  # colnames(E) <- paste0("K_", 1:5)
  # labels = factor(rep(letters[1:2], each = 10))

  if (is.null(labels) & is.null(local)) {
    labels = factor(rep("All", nrow(E)))
  }

  if (!is.null(labels)) {
    if (class(labels) != "factor") {
      labels <- factor(labels)
    }
    L = Matrix::fac2sparse(labels)

    LE = L %*% E

    k_best = getArgMin(LE, ...)

    if (outputPerCell) {
      k_best <- k_best[labels]
      names(k_best) <- rownames(E)
    }

    return(k_best)
  }

  # if function still running, then use the neighbours in local
  # ensure that self is also included
  # local_self = cbind(seq_len(nrow(E)), local)
  local_self = local

  LE = apply(E, 2, function(e) Matrix::rowSums(vectorSubset(e, local_self)))

  k_best = getArgMin(LE, ...)
  names(k_best) <- rownames(E)

  return(k_best)
}
