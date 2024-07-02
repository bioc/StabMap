test_that("stabMap works", {
  set.seed(2021)
  assay_list <- mockMosaicData()

  # specify which datasets to use as reference coordinates
  reference_list <- c("D1", "D3")

  # specify some sample labels to distinguish using linear discriminant
  # analysis (LDA)
  labels_list <- list(
    D1 = rep(letters[1:5], length.out = ncol(assay_list[["D1"]]))
  )

  # stabMap
  out <- stabMap(assay_list,
    reference_list = reference_list,
    labels_list = labels_list,
    ncomponentsReference = 20,
    ncomponentsSubset = 20,
    verbose = FALSE
  )

  testthat::expect_snapshot_value(
    out,
    style = "serialize"
  )
})
