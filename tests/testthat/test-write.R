test_that("write csv works as expected", {
  fx_dir <- system.file("extdata", "hps40", package = "hirisplexr")
  prefix <- file.path(fx_dir, "hps40")

  exp1 <- c(rep(2L, 10), rep(1L, 10), rep(0L, 10), rep(2L, 10), NA)

  out1 <- tempfile(fileext = ".csv")
  write_hirisplex_csv(prefix, panel = "hirisplexs", out = out1,
                      sample_id = "IID", allow_strand_flip = TRUE)
  df <- read.csv(out1,na.strings="NA")
  row1 <- as.integer(as.vector(unname(df[1,2:42])))
  
  testthat::expect_equal(row1, exp1)
})
