test_that("panel CSV loads and has expected columns", {
  loader <- getFromNamespace(".load_hirisplex_panels", "hirisplexr")
  p <- loader()
  expect_true(all(c("panel","order","rsid","input_allele","column_name") %in% names(p)))
  expect_true(any(p$panel == "hirisplexs"))
  expect_true(any(p$panel == "hirisplex"))
  expect_true(any(p$panel == "irisplex"))
  # Ordering must start at 1 for each panel
  expect_equal(min(p$order[p$panel=="hirisplexs"]), 1)
})
