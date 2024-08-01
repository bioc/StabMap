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

  # error when out is not provided
  testthat::expect_error(quietIE(assay_list))

  # error when rownames or colnames not given
  out2 <- out
  rownames(out2) <- NULL
  colnames(out2) <- NULL
  testthat::expect_error(quietIE(assay_list, out2))

  # expect same output when out is sparse
  out_sparse = as(out, "dgCMatrix")
  imp_sparse = imp <- quietIE(assay_list, out_sparse)
  testthat::expect_equal(imp, imp_sparse)
})
