#' Adaptive k-Nearest Neighbour Classification
#'
#' Adaptive k-Nearest Neighbour Classification for a k-nearest neighbour matrix,
#' given class labels and local k values for the training data
#'
#'
#' @param knn Is a k-nearest neighbour matrix, giving the indices of the
#' training set that the query is closest to. Rows are the query cells, columns
#' are the NNs. Typically output using
#' BiocNeighbors::queryKNN(,,k = max(k_local)).
#' @param class Is the labels associated with the training set.
#' @param k_local Is an integer vector length of the training set, giving the
#' local k to use if k_local is given as a single integer, then that value is
#' used as k for all observations.
#'
#' @return A character vector of of classifications for the test set.
#'
#' @examples
#' # Generate example data
#' # data = matrix(rpois(10*20, 10), 10, 20) # 10 genes, 20 cells
#' # data_2 = matrix(rpois(10*30, 10), 10, 30) # 10 genes, 30 cells
#'
#' # Generate error matrix for k_local
#' # E = matrix(runif(100),20,5)
#' # colnames(E) <- paste0("K_", 1:5)
#'
#' # Define training class labels and adaptive k-values
#' # class = factor(rep(letters[1:2], each = 10))
#' # k_local = getAdaptiveK(E, labels = class)
#'
#' # knn = queryKNN(t(data), t(data_2), k = max(as.numeric(gsub("K_", "", k_local))))$index
#'
#' # Adaptive KNN classification
#' # test = adaptiveKNN(knn, class, as.numeric(gsub("K_", "", k_local)))
#'
#' @export
adaptiveKNN = function(knn,
                       class,
                       k_local) {

  # knn is a k-nearest neighbour matrix, giving the
  # indices of the training set that the query is
  # closest to. Rows are the query cells, columns
  # are the NNs, should be a large value. Typically
  # output using BiocNeighbors::queryKNN(,,k = max(k_local))

  # class is the labels associated with the training
  # set

  # k_local is an integer vector length of the training
  # set, giving the local k to use
  # if k_local is given as a single integer, then
  # that value is used as k for all observations

  # example data
  # data = matrix(rpois(10*20, 10), 10, 20) # 10 genes, 20 cells
  # local = BiocNeighbors::findKNN(t(data), k = 5, get.distance = FALSE)$index
  # A = matrix(runif(100),20,5)
  # colnames(A) <- paste0("K_", 1:5)
  # labels = factor(rep(letters[1:2], each = 10))
  # k_local = getAdaptiveK(A, labels = labels)
  # data_2 = matrix(rpois(10*30, 10), 10, 30) # 10 genes, 30 cells
  # knn = BiocNeighbors::queryKNN(t(data), t(data_2), k = 5, get.distance = FALSE)$index
  # class = labels

  # if (length(k_local) == 1) {
  #   k_local <- rep(k_local, max(knn))
  # }
  #
  # # first use 1NN to identify the local best k value
  # query_best_k = k_local[knn[,1]]

  query_best_k = getQueryK(knn, k_local)

  if (any(query_best_k > ncol(knn))) {
    warning("k is larger than nearest neighbours provided, taking all neighbours given")
  }

  # convert the KNN names to the class labels
  knn_class = vectorSubset(class, knn)

  # extract the most frequent among the nearest local best k neighbours
  # new_class = sapply(1:nrow(knn), function(i) getModeFirst(knn_class[i,], query_best_k[i]))
  # same as:
  new_class = mapply(getModeFirst, split(knn_class, seq_len(nrow(knn_class))), query_best_k,
                     USE.NAMES = FALSE)
  names(new_class) <- rownames(knn_class)

  return(new_class)
}
