# Internal: load the packaged panel definitions
# Stored as inst/extdata/hirisplex_panels.csv
.load_hirisplex_panels <- function() {
  path <- system.file("extdata", "hirisplex_panels.csv", package = "hirisplexr", mustWork = TRUE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  # Basic checks
  stopifnot(all(c("panel","order","rsid","input_allele","column_name") %in% names(x)))
  # Ensure ordering by 'order' within each panel
  x <- x[order(x$panel, x$order), ]
  x
}
