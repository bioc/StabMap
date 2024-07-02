#' getArgMin
#'
#' For each row in a matrix calculate the first index which gives the minimum
#' value
#'
#' @param M A matrix.
#' @param return_colnames Logical whether to return column names of matrix
#' (default TRUE). Otherwise return index.
#' @param identicalNA Logical whether to return NA if all values in a row are
#' identical (default TRUE).
#'
#' @return A vector containing the first index or column name of the
#' minimum values for each row of the matrix.
#'
#' @keywords internal
getArgMin <- function(M, return_colnames = TRUE, identicalNA = TRUE) {
  m <- max.col(-M, ties.method = "first")

  if (return_colnames) {
    if (!is.null(colnames(M)[1])) {
      m <- colnames(M)[m]
    }
  }

  names(m) <- rownames(M)

  if (identicalNA) {
    # if all the values in a row of M are identical,
    # return NA
    m[apply(M, 1, allEqual)] <- NA
  }

  return(m)
}


#' getModeFirst
#'
#' Identify the mode of x up to the first index
#'
#' @param x A character or a factor.
#' @param first An integer.
#'
#' @return A character of the mode of x.
#'
#' @keywords internal
getModeFirst <- function(x, first) {
  names(which.max(table(x[seq_len(first)])[unique(x[seq_len(first)])]))
}

#' getQueryK
#'
#' For each cell in the query data, use the 1NN's adaptive k value
#' (of the reference data) to identify the local best k value
#'
#' @param knn Is a k-nearest neighbour matrix, giving the indices of the
#' training set that the query is closest to. Rows are the query cells, columns
#' are the NNs, should be a large value. Typically output using
#' BiocNeighbors::queryKNN(,,k = max(k_local)).
#' @param k_local Is an integer vector length of the reference set, giving the
#' local k to use. If k_local is given as a single integer, then that value is
#' used as k for all observations.
#'
#' @return An integer vector with local k to use for each query cell.
#'
#' @keywords internal
getQueryK <- function(knn, k_local) {

  if (length(k_local) == 1) {
    k_local <- rep(k_local, nrow(knn))
    return(k_local)
  }

  query_best_k <- k_local[knn[, 1]]

  return(query_best_k)
}


#' isUnequal
#'
#' Checks if elements of 2 vectors are unequal
#'
#' @param x A vector.
#' @param y A vector.
#'
#' @return An integer vector. 1 for unequal. 0 for equal
#'
#' @keywords internal
isUnequal <- function(x, y) {
  1 * (as.character(x) != as.character(y))
}

#' allEqual
#'
#' Checks if a vector is equal to its first element
#'
#' @param x A vector.
#'
#' @return logical whether a a vector is equal to its first element.
#'
#' @keywords internal
allEqual <- function(x) {
  all(x == x[1])
}

#' getBinaryError
#'
#' For potential k values, generate a binary error matrix from KNN label
#' classification
#'
#' @param knn Is a k-nearest neighbour matrix, giving the indices of the
#' training set that the query is closest to. Rows are the query cells, columns
#' are the NNs, should be a large value. Typically output using
#' BiocNeighbors::queryKNN(,,k = max(k_values)).
#' @param k_values Is an integer vector of the values of k to consider for
#' extracting accuracy. If k_values has names then pass these to colnames of E.
#' @param class_train Is a factor or character vector of classes that
#' corresponds to the indices given within knn.
#' @param class_true Is a factor or character vector that corresponds to the
#' rows of knn. If class_true has names then pass these to rownames of E.
#'
#' @return A sparse binary error matrix.
#'
#' @keywords internal
getBinaryError <- function(knn,
                           k_values,
                           class_train,
                           class_true) {

  if (max(unlist(k_values)) > ncol(knn)) {
    m <- paste0(
      "largest k value is larger than neighbours provided in knn (i.e. ",
      ncol(knn),
      "), select a smaller k value or provide more neighbours"
    )
    stop(m)
  }

  class_pred <- mapply(
    adaptiveKNN, k_values,
    MoreArgs = list(class = class_train, knn = knn)
  )
  E <- methods::as(
    apply(class_pred, 2, isUnequal, y = class_true), "sparseMatrix"
  )

  rownames(E) <- names(class_true)
  colnames(E) <- names(k_values)

  return(E)
}

#' combineBinaryErrors
#'
#' Combines binary error matrices by averaging error values across all matrices,
#' for each entry (row and column combination)
#'
#' @param E_list A list containing matrices. Each matrix must have the same
#' number of columns (k-values) and contain rownames (cells).
#'
#' @return A sparse error matrix.
#'
#' @keywords internal
combineBinaryErrors <- function(E_list) {

  if (!allEqual(unlist(lapply(E_list, ncol)))) {
    stop("each matrix in E_list must have the same number of columns")
  }

  if (any(unlist(lapply(E_list, function(x) is.null(rownames(x)))))) {
    stop("each matrix in E_list must have associated rownames")
  }

  all_rows <- unlist(lapply(E_list, rownames))

  E_exp <- do.call(rbind, E_list)
  E_split <- split.data.frame(E_exp, all_rows)

  # data becomes dense at this stage
  E_means <- lapply(E_split, Matrix::colMeans, na.rm = TRUE)

  E <- methods::as(do.call(rbind, E_means), "sparseMatrix")

  return(E)
}


#' getBestColumn
#'
#' Identifies the index of the column of a matrix with the minimum mean.
#' If balanced_labels is given then calculate the balanced mean
#'
#' @param E An error matrix.
#' @param balanced_labels Class labels for each row (cell) of E.
#'
#' @return The index of the best performing column of E
#'
#' @keywords internal
getBestColumn <- function(E,
                          balanced_labels = NULL) {

  if (is.null(balanced_labels)) {
    return(which.min(colMeans(E, na.rm = TRUE)))
  } else {
    E_lab <- apply(E, 2, function(x) tapply(x, balanced_labels, mean))
    return(which.min(colMeans(E_lab, na.rm = TRUE)))
  }
}

#' smoothLocal
#'
#' Smooth out the adaptive k values. Can be smoothed by computing the
#' arithmetic or geometric mean of the adaptive k-values for each cells
#' neighbourhood
#'
#' @param best_k Is a named vector of local best k values
#' @param local Is a KNN matrix, with rows same as best_k and values indices of
#' best_k.
#' @param smooth An integer of k-nearest neighbours to smooth over.
#' @param mean_type Character indicating to calculate the 'geometric' or
#' 'arithmetic' mean.
#'
#' @return A numeric vector of smoothed adaptive k-values.
#'
#' @keywords internal
smoothLocal <- function(best_k, local, smooth = 10, mean_type = "geometric") {

  nms <- names(best_k)

  if (mean_type == "arithmetic") {
    best_k_smooth <- ceiling(
      rowMeans(vectorSubset(best_k, local[, seq_len(smooth)]))
    )
  }
  if (mean_type == "geometric") {
    best_k_smooth <- ceiling(
      apply(vectorSubset(best_k, local[, seq_len(smooth)]), 1, gm_mean)
    )
  }
  names(best_k_smooth) <- nms
  return(best_k_smooth)
}

#' gm_mean
#'
#' Calculate the geometric mean
#'
#' @param x A vector.
#' @param na.rm A logical value indicating whether NA values should be stripped
#' before calculating the geometric mean.
#'
#' @return A numeric.
#'
#' @keywords internal
gm_mean <- function(x, na.rm = TRUE) {
  if (all(is.na(x))) {
    return(NA)
  }
  exp(sum(log(x[x > 0]), na.rm = na.rm) / sum(!is.na(x)))
}

#' buildLabelsDataFrame
#'
#' Build dataframe for output from `classifyEmbedding`
#'
#' @param labels Is a named character vector with true labels.
#' @param resubstituted_labels Is a named character vector with predicted
#' labels.
#' @param k_adaptive Is a named vector of the k-values, this could be a single
#' integer when fixed.
#'
#' @return A dataframe with rows the same as resubstituted_labels and columns
#' for input_labels, predicted_labels, and resubstituted_labels.
#'
#' @keywords internal
buildLabelsDataFrame <- function(labels, resubstituted_labels, k_adaptive) {

  df <- data.frame(
    input_labels = as.character(labels[names(resubstituted_labels)]),
    resubstituted_labels = as.character(resubstituted_labels),
    predicted_labels = as.character(resubstituted_labels),
    row.names = names(resubstituted_labels)
  )
  df[names(labels), "predicted_labels"] <- labels

  if (length(k_adaptive) == 1) {
    df$k <- k_adaptive
  } else {
    df$k <- k_adaptive[rownames(df)]
  }

  return(df)
}

#' getBinaryErrorFromPredictions
#'
#' Compute binary error between predicted labels and true labels
#'
#' @param pred Is a matrix of class label predictions.
#' @param labels Is a named vector of true labels.
#'
#' @return A sparse binary error matrix.
#'
#' @keywords internal
getBinaryErrorFromPredictions <- function(pred, labels) {

  E <- methods::as(
    apply(pred, 2, isUnequal, y = labels[rownames(pred)]), "sparseMatrix"
  )
  rownames(E) <- names(labels)
  colnames(E) <- colnames(pred)

  return(E)
}

#' queryNamedKNN
#'
#' queryNamedKNN
#'
#' @param coords_reference coords_reference
#' @param coords_query coords_query
#' @param k k
#'
#' @return matrix
#'
#' @keywords internal
queryNamedKNN <- function(coords_reference, coords_query, k) {

  knn <- BiocNeighbors::queryKNN(
    coords_reference,
    coords_query,
    k = k, get.distance = FALSE
  )$index
  rownames(knn) <- rownames(coords_query)
  knn_name <- vectorSubset(rownames(coords_reference), knn)

  return(knn_name)
}
