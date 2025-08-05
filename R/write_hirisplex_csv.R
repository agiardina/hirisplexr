#' Write HIrisPlex / HIrisPlex-S CSV from a PLINK BED/BIM/FAM prefix
#'
#' Given a PLINK 1.9 binary dataset (prefix.bed/.bim/.fam), this function
#' produces a CSV ready to upload to the HIrisPlex(-S) webtool.
#'
#' Columns are 'SampleID' followed by one column per required SNP in the form
#' `rsID_Allele` (e.g., `rs12203592_T`). Each cell contains 0/1/2 (count of the
#' input allele) or NA when the SNP is missing. The column set and order are
#' defined by the selected panel (IrisPlex, HIrisPlex, HIrisPlex-S).
#'
#' Allele counting is based on the PLINK .bim alleles. Genotype dosage is read
#' on demand from the .bed using [BEDMatrix::BEDMatrix()], which encodes the
#' dosage of the first allele in the .bim file (A1). If the webtool's required
#' input allele equals A1, we use the dosage directly; if it equals A2, we use
#' (2 - dosage). If `allow_strand_flip = TRUE`, we also reconcile complements
#' (A<->T, C<->G) to account for strand orientation differences.
#'
#' @param prefix Character. Path prefix to PLINK files, without extension.
#' @param panel Character. One of "hirisplexs" (default), "hirisplex", "irisplex".
#' @param out Character. Output CSV path. Defaults to `<panel>.csv` in the
#'   working directory.
#' @param sample_id Character. How to form 'SampleID': "IID" or "FID_IID".
#'   Default is "IID".
#' @param allow_strand_flip Logical. If TRUE, attempt to match the required
#'   allele by allowing strand complements.
#' @return (Invisibly) the output file path.
#' @examples
#' \dontrun{
#' write_hirisplex_csv("/path/to/prefix", panel = "hirisplexs")
#' }
#' @export
write_hirisplex_csv <- function(prefix,
                                panel = c("hirisplexs", "hirisplex", "irisplex"),
                                out = NULL,
                                sample_id = c("IID", "FID_IID"),
                                allow_strand_flip = TRUE) {
  panel     <- match.arg(panel)
  sample_id <- match.arg(sample_id)

  if (is.null(out)) out <- paste0(panel, ".csv")

  # --- Load panel specification (rsid, allele to count, column_name) ---
  p <- .load_hirisplex_panels()
  p <- p[p$panel == panel, , drop = FALSE]

  if (nrow(p) == 0L) stop("Unknown panel: ", panel)

  # --- Read BIM/FAM and map indices for required rsIDs ---
  bim_path <- paste0(prefix, ".bim")
  fam_path <- paste0(prefix, ".fam")
  bed_path <- paste0(prefix, ".bed")

  if (!file.exists(bim_path) || !file.exists(fam_path) || !file.exists(bed_path)) {
    stop("Missing PLINK files. Expected ", shQuote(bed_path), ", ",
         shQuote(bim_path), " and ", shQuote(fam_path), ".")
  }

  bim <- data.table::fread(bim_path,
    col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  fam <- data.table::fread(fam_path,
    col.names = c("FID", "IID", "PID", "MID", "SEX", "PHENO"))

  # Vector of rsIDs in BIM and map from rsID -> index
  idx_in_bim <- match(p$rsid, bim$SNP)  # preserves panel order
  
  # Memory-map genotypes
  G <- BEDMatrix::BEDMatrix(bed_path)

  # Build SampleID column
  SampleID <- if (sample_id == "IID") fam$IID else paste(fam$FID, fam$IID, sep = "_")

  # Prepare output matrix (integer, with NAs)
  out_mat <- matrix(NA_integer_, nrow = nrow(fam), ncol = nrow(p))
  colnames(out_mat) <- p$column_name

  # Helper: base complement
  complement <- function(a) chartr("ACGT", "TGCA", a)

  # For each panel SNP, compute allele counts (0/1/2) or NA if missing/unreconcilable
  for (j in seq_len(nrow(p))) {
    i <- idx_in_bim[j]
    need <- p$input_allele[j]

    if (is.na(i)) {
      # SNP missing in BIM: leave NA
      next
    }

    a1 <- bim$A1[i]; a2 <- bim$A2[i]
    
    # Dosage of A1 from BEDMatrix, as integer 0/1/2 with NAs
    dosA1 <- G[, i, drop = FALSE]

    use_direct <- identical(need, a1) ||
      (allow_strand_flip && identical(complement(need), a1))

    use_flip2 <- identical(need, a2) ||
      (allow_strand_flip && identical(complement(need), a2))

    if (isTRUE(use_direct)) {
      out_mat[, j] <- dosA1
    } else if (isTRUE(use_flip2)) {
      out_mat[, j] <- 2L - dosA1
    } else {
      # Unreconcilable allele labels -> leave NA, issue a warning
      warning(sprintf("Alleles not reconcilable for %s (A1=%s, A2=%s, need=%s)",
                      p$rsid[j], a1, a2, need))
    }
  }

  # Bind and write
  out_dt <- data.table::data.table(SampleID = SampleID, out_mat)
  data.table::fwrite(out_dt, out, na = "NA")

  invisible(out)
}
