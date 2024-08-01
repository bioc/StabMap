#' Adaptive k-Nearest Neighbour Classification using the StabMap joint embedding
#'
#' Performs adaptive k-nearest neighbour classification of discrete labels for a
#' training set from a query set, leveraging the StabMap joint embedding. The
#' training labels are defined in `labels`, with all other rows of the
#' embedding treated as the testing set.
#'
#' @param coords A cells (rows) x dimensions data matrix, on which euclidean
#' distances are to be calculated for KNN classification. Must have rownames.
#' Typically, output from `stabMap()`.
#' @param labels A named character vector of labels for the training set.
#' @param type A character of the type of adaptive KNN classification to be
#' used. Must be one of "adaptive_local", "adaptive_labels",
#' "uniform_optimised", or "uniform_fixed". Default
#' is "uniform_fixed".
#' @param k_values A numeric vector of potential k values. If type is
#' "uniform_fixed", then the first value of k_values is used. Default is 5.
#' @param error_measure Is the error type to use for selection of the best k.
#' Must be one of "simple_error" or "balanced_error". "simple_error" weights all
#' cells equally. "balanced_error" weights error by `labels` factors. Only
#' affects error type for type == "uniform_optimised".
#' @param adaptive_nFold Is the number of folds for adaptive selection
#' cross-validation.
#' @param adaptive_nRep Is the number of repetitions of adaptive selection
#' cross-validation.
#' @param adaptive_local_nhood Is the neighbourhood size for optimising locally.
#' @param adaptive_local_smooth Is the number of neighbours to use for smoothing
#' locally.
#' @param verbose Logical whether to print repetition and fold number for
#' adaptive selection cross-validation.
#'
#' @return Is a dataframe with rows the same as coords, and same rownames.
#' Columns are: input_labels is the training labels that were provided in
#' `labels` (NA is used as labels for the testing set), resubstituted_labels is
#' predicted labels for all rows (including for the training data),
#' predicted_labels is predicted labels for the testing set but true labels as
#' provided in `labels` for the training set, k is the adaptive k value used for
#' that each row of the training set.
#'
#' @examples
#' set.seed(100)
#' # Simulate coordinates
#' coords <- matrix(rnorm(1000), 100, 10)
#' rownames(coords) <- paste0("cell_", 1:nrow(coords))
#'
#' # Define labels of the first 50 cells
#' labels <- rep(paste0("type_", letters[1:5]), 10)
#' names(labels) <- rownames(coords)[1:length(labels)]
#'
#' # Uniform fixed KNN classification
#' knn_out <- classifyEmbedding(
#'   coords, labels,
#'   type = "uniform_fixed", k_values = 5
#' )
#' table(knn_out$predicted_labels)
#'
#' # Adaptive KNN classification using local error
#' knn_out <- classifyEmbedding(
#' coords, labels,
#' type = "adaptive_local",
#' k_values = 2:3,
#' adaptive_nFold = 5,
#' adaptive_nRep = 10
#' )
#' table(knn_out$predicted_labels)
#'
#  # Adaptive KNN classification using adaptive labels
#' knn_out <- classifyEmbedding(
#'   coords, labels,
#'   type = "adaptive_labels",
#'   k_values = 2:3,
#'   adaptive_nFold = 5,
#'   adaptive_nRep = 10
#' )
#' table(knn_out$predicted_labels)
#'
#' # Adaptive KNN classification using uniform optimised with balanced error
#' knn_out <- classifyEmbedding(
#'   coords, labels,
#'   type = "uniform_optimised",
#'   k_values = 2:3,
#'   adaptive_nFold = 3,
#'   adaptive_nRep = 10,
#'   error_measure = "balanced_error"
#' )
#' table(knn_out$predicted_labels)
#'
#' @export
classifyEmbedding <- function(
    coords,
    labels,
    type = c("uniform_fixed", "adaptive_labels", "adaptive_local", "uniform_optimised"),
    k_values = 5,
    error_measure = c("simple_error", "balanced_error"),
    adaptive_nFold = 2,
    adaptive_nRep = 5,
    adaptive_local_nhood = 100,
    adaptive_local_smooth = 10,
    verbose = TRUE) {

  type <- match.arg(type)
  error_measure <- match.arg(error_measure)

  if (is.null(rownames(coords)[1])) {
    stop("coords must have rownames")
  }

  if (is.null(names(labels)[1])) {
    stop("labels must have names")
  }

  max_k <- max(k_values)

  coords_train <- coords[names(labels), ]

  if (type == "uniform_fixed") {
    k <- k_values[1]

    knn <- queryNamedKNN(coords_train, coords, k = k)

    resubstituted_labels <- adaptiveKNN(
      knn,
      labels,
      k
    )

    out <- buildLabelsDataFrame(labels, resubstituted_labels, k_adaptive = k)

    return(out)
  }

  E_list <- vector("list", length = adaptive_nRep)

  for (Rep in seq_len(adaptive_nRep)) {
    train_k <- sample(seq_len(adaptive_nFold), nrow(coords_train),
      prob = rep(1, adaptive_nFold), replace = TRUE
    )

    if (verbose) message("Rep ", Rep, " of ", adaptive_nRep)

    for (fold in seq_len(adaptive_nFold)) {
      if (verbose) message("Fold ", fold, " of ", adaptive_nFold)

      coords_train_k <- coords_train[train_k == fold, ]
      coords_test_k <- coords_train[!train_k == fold, ]

      knn <- queryNamedKNN(coords_train_k, coords_test_k, k = max_k)

      labels_train <- labels[rownames(coords_train_k)]
      labels_test <- labels[rownames(coords_test_k)]

      E <- getBinaryError(
        knn = knn,
        k_values = k_values,
        class_train = labels_train,
        class_true = labels_test
      )

      E_list[[Rep]] <- E
    }
  }
  E <- combineBinaryErrors(E_list)[rownames(coords_train), ]


  if (type == "uniform_optimised") {
    # select the best column based on the error type
    if (error_measure == "simple_error") {
      best_column <- getBestColumn(E)
    }
    if (error_measure == "balanced_error") {
      best_column <- getBestColumn(E, balanced_labels = labels[rownames(E)])
    }
    k_adaptive <- k_values[best_column]
  }

  if (type == "adaptive_labels") {
    best_k_labels_index <- getAdaptiveK(E,
      labels = labels[rownames(E)],
      local = NULL,
      outputPerCell = TRUE
    )
    best_k_labels <- vectorSubset(k_values, as.matrix(best_k_labels_index))[, 1]

    # if any are NA, select the geometric mean value among the rest
    best_k_labels[is.na(best_k_labels)] <- ceiling(gm_mean(best_k_labels))

    k_adaptive <- best_k_labels[names(labels)]
  }

  if (type == "adaptive_local") {
    local <- queryNamedKNN(coords_train, coords_train, k = adaptive_local_nhood)

    best_k_local_index <- getAdaptiveK(E,
      local = local,
      outputPerCell = TRUE
    )
    best_k_local_unsmoothed <- vectorSubset(
      k_values, as.matrix(best_k_local_index)
    )[, 1]

    if (any(is.na(best_k_local_unsmoothed))) {
      defined <- names(best_k_local_unsmoothed)[!is.na(best_k_local_unsmoothed)]

      local_defined <- queryNamedKNN(
        coords_train[defined, ], coords_train,
        k = adaptive_local_nhood
      )
    } else {
      local_defined <- local
    }

    best_k_local <- smoothLocal(
      best_k_local_unsmoothed, local_defined,
      smooth = adaptive_local_smooth
    )

    k_adaptive <- best_k_local[names(labels)]
  }

  max_k <- max(k_adaptive)

  knn <- queryNamedKNN(coords_train, coords, k = max_k)

  resubstituted_labels <- adaptiveKNN(
    knn,
    labels,
    k_adaptive
  )

  out <- buildLabelsDataFrame(
    labels, resubstituted_labels,
    k_adaptive = k_adaptive
  )

  return(out)
}
