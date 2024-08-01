test_that("classifyEmbedding works", {
  set.seed(1234)

  coords <- matrix(rnorm(1000), 100, 10)
  rownames(coords) <- paste0("cell_", seq_len(nrow(coords)))

  # Define labels of the first 50 cells
  labels <- rep(paste0("type_", letters[1:5]), 10)
  names(labels) <- rownames(coords)[seq_along(labels)]

  quietCE <- purrr::quietly(classifyEmbedding)

  # Uniform fixed KNN classification
  knn_out <- quietCE(
    coords, labels,
    type = "uniform_fixed", k_values = 5
  )
  testthat::expect_snapshot_value(
    knn_out, style = "serialize"
  )

  # Adaptive KNN classification using local error
  knn_out <- quietCE(
    coords, labels,
    type = "adaptive_local",
    k_values = 2:3,
    adaptive_nFold = 5,
    adaptive_nRep = 10
  )
  testthat::expect_snapshot_value(
    knn_out, style = "serialize"
  )

  # Adaptive KNN classification using adaptive labels
  knn_out <- quietCE(
    coords, labels,
    type = "adaptive_labels",
    k_values = 2:3,
    adaptive_nFold = 5,
    adaptive_nRep = 10
  )

  testthat::expect_snapshot_value(
    knn_out, style = "serialize"
  )

  # Adaptive KNN classification using uniform optimised with balanced error
  knn_out <- quietCE(
    coords, labels,
    type = "uniform_optimised",
    k_values = 2:3,
    adaptive_nFold = 3,
    adaptive_nRep = 10,
    error_measure = "balanced_error"
  )

  testthat::expect_snapshot_value(
    knn_out, style = "serialize"
  )
})
