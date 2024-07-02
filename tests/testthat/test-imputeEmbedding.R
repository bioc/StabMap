test_that("imputeEmbedding works", {
  set.seed(2021)
  assay_list <- mockMosaicData()

  # stabMap
  out <- stabMap(assay_list,
    ncomponentsReference = 20,
    ncomponentsSubset = 20,
    verbose = FALSE
  )

  # impute values
  quietIE <- purrr::quietly(imputeEmbedding)
  imp <- quietIE(assay_list, out)

  testthat::expect_snapshot_value(
    imp, style = "serialize"
  )

})
