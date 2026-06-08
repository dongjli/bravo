## ----------------------------------------------------------------------
## Auto-detect the field separator from the first few lines of a file.
## Returns list(sep = <string>, fixed = <logical>).
## ----------------------------------------------------------------------
detect_separator <- function(file.name, n_peek = 10L) {
  con <- file(file.name, open = "r")
  on.exit(close(con))
  lines <- readLines(con, n = n_peek)
  lines <- lines[nzchar(lines)]
  if (length(lines) < 1L) stop("File appears to be empty: ", file.name)
  
  ## pattern + whether it is a fixed string or a regular expression
  candidates <- list(
    list(sep = "\t",   fixed = TRUE),    # tab
    list(sep = ",",    fixed = TRUE),    # comma
    list(sep = ";",    fixed = TRUE),    # semicolon
    list(sep = " ",    fixed = TRUE),    # single space
    list(sep = "\\s+", fixed = FALSE)    # runs of whitespace (fallback)
  )
  
  score <- vapply(candidates, function(cd) {
    counts <- vapply(
      lines,
      function(l) length(strsplit(trimws(l), cd$sep, fixed = cd$fixed)[[1]]),
      integer(1)
    )
    ## a good separator splits every peeked line into the SAME (>1) number
    ## of fields; score = that field count, or 0 if inconsistent / absent
    if (counts[1] > 1L && all(counts == counts[1])) counts[1] else 0L
  }, integer(1))
  
  if (max(score) <= 1L)
    stop("Could not auto-detect a separator; please pass `separator=` explicitly.")
  
  best <- candidates[[which.max(score)]]
  pretty <- switch(best$sep, "\t" = "<tab>", " " = "<space>",
                   "\\s+" = "<whitespace>", best$sep)
  message("Auto-detected separator: ", pretty)
  best
}

## ----------------------------------------------------------------------
## Convert a dense numeric genotype file to a sparse Matrix.
##
##   file.name     path to the file (header row, then "name <vals...>" rows)
##   num.genotypes optional; if NULL it is counted in pass 1
##   separator     optional; auto-detected if NULL
##   fixed         is `separator` a fixed string (TRUE) or a regex (FALSE)?
##   recode        recode each SNP so its major allele is 0 (per-SNP)
##   na.strings    tokens treated as missing
##   repr          "C" (dgCMatrix, default) or "R" (dgRMatrix)
## ----------------------------------------------------------------------

#' @title Convert numeric genotype matrix to sparse matrix
#' @description Reads a numeric genotype file and converts it to a sparse matrix format.
#' @param file.name Path to the numeric genotype file. Could be (and should be) gzipped.
#' @param num.genotypes Maximum number of genotypes to read. An upper bound is OK.
#' @param separator "\\t" or "," etc. that separates the entries in a line.
#' @param fixed Is `separator` a fixed string (TRUE) or a regex (FALSE)?
#' @param recode        recode each SNP so its major allele is 0 (per-SNP)
#' @param na.strings    tokens treated as missing
#' @param progress Whether to show a progress bar (default TRUE).
#' @return A sparse matrix of class dgCMatix.
#' @export
dense2sparse <- function(file.name,
                         num.genotypes = NULL,
                         separator     = NULL,
                         fixed         = TRUE,
                         recode        = TRUE,
                         na.strings    = c("NA", ".", "-9", "NN"),
                         progress      = TRUE) {
  
  # repr <- match.arg(repr)
  repr <- "C" # In the future we might support "R".
  
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("The 'Matrix' package is required.")
  
  ## --- resolve the separator -------------------------------------------
  if (is.null(separator)) {
    sp        <- detect_separator(file.name)
    separator <- sp$sep
    fixed     <- sp$fixed
  }
  split_line <- function(s) strsplit(trimws(s), separator, fixed = fixed)[[1]]
  to_num <- function(ch) {
    ch[ch %in% na.strings] <- NA
    suppressWarnings(as.numeric(ch))
  }
  
  ## --- read the header -------------------------------------------------
  con    <- file(file.name, open = "r")
  header <- readLines(con, n = 1L)
  close(con)
  col_names <- split_line(header)[-1L]
  n_snp     <- length(col_names)
  if (n_snp < 1L) stop("No SNP columns found in the header.")
  
  ## =====================================================================
  ## PASS 1: count genotypes and, if recoding, accumulate per-SNP means.
  ## Skipped only when recoding is off AND the row count is already known.
  ## =====================================================================
  flip   <- logical(n_snp)
  n_geno <- num.genotypes
  
  ploidy = -1;
  
  if (recode || is.null(num.genotypes)) {
    col_sum <- numeric(n_snp)
    col_n   <- numeric(n_snp)
    n_geno  <- 0L
    
    con <- file(file.name, open = "r")
    invisible(readLines(con, n = 1L))                       # skip header
    repeat {
      line <- readLines(con, n = 1L)
      if (length(line) == 0L) break
      if (!nzchar(trimws(line))) next
      if (!is.null(num.genotypes) && n_geno >= num.genotypes) break
      
      vals <- to_num(split_line(line)[-1L])
      if (length(vals) != n_snp)
        stop(sprintf("Data row %d has %d values but the header has %d SNPs.",
                     n_geno + 1L, length(vals), n_snp))
      n_geno <- n_geno + 1L
      
      if(any(is.na(vals)))
        stop("Missing values (NA) found. Please impute missing values.")

      if (recode) {
        ploidy <- max(ploidy, max(vals,na.rm = TRUE))
        col_sum <- col_sum + vals
        col_n   <- col_n   + 1
      }
    }
    close(con)
    if (n_geno == 0L) stop("No genotype rows found.")
    
    ## flip a SNP when its mean dosage exceeds ploidy/2 -> major allele -> 0
    if (recode) {
      col_mean <- ifelse(col_n > 0, col_sum / col_n, 0)
      flip     <- col_mean > (ploidy / 2)
      if (any(flip))
        message(sprintf("Recoding %d of %d SNP(s) so the major allele is 0.",
                        sum(flip), n_snp))
    }
  }
  
  ## =====================================================================
  ## PASS 2: build the matrix, storing only non-zero (and NA) entries.
  ## =====================================================================
  row_names   <- character(n_geno)
  col_indices <- vector("list", n_geno)   # 0-based column indices per row
  row_data    <- vector("list", n_geno)
  row_ptr     <- integer(n_geno + 1L)
  
  if (progress) pb <- txtProgressBar(min = 0, max = n_geno, style = 3)
  
  con <- file(file.name, open = "r")
  invisible(readLines(con, n = 1L))                         # skip header
  i <- 0L
  repeat {
    line <- readLines(con, n = 1L)
    if (length(line) == 0L || i >= n_geno) break
    if (!nzchar(trimws(line))) next
    
    i      <- i + 1L
    fields <- split_line(line)
    row_names[i] <- fields[1L]
    vals   <- to_num(fields[-1L])
    if (length(vals) != n_snp)
      stop(sprintf("Data row %d has %d values but the header has %d SNPs.",
                   i, length(vals), n_snp))
    if (any(flip)) vals[flip] <- ploidy - vals[flip]
    
    ## keep non-zeros AND NAs so missing data is never silently zero
    nz <- which(vals != 0 | is.na(vals))
    col_indices[[i]] <- nz - 1L
    row_data[[i]]    <- vals[nz]
    row_ptr[i + 1L]  <- row_ptr[i] + length(nz)
    if (progress) setTxtProgressBar(pb, i)
  }
  close(con)
  if (progress) close(pb)
  
  ## trim if the file held fewer rows than expected
  if (i < n_geno) {
    n_geno      <- i
    row_names   <- row_names[seq_len(i)]
    col_indices <- col_indices[seq_len(i)]
    row_data    <- row_data[seq_len(i)]
    row_ptr     <- row_ptr[seq_len(i + 1L)]
  }
  message("\nGenotypes read: ", n_geno)
  
  ## --- assemble (j + p => row-compressed input; coerced to `repr`) -----
  M <- Matrix::sparseMatrix(
    j        = unlist(col_indices, use.names = FALSE),
    p        = row_ptr,
    x        = unlist(row_data,    use.names = FALSE),
    dims     = c(n_geno, n_snp),
    dimnames = list(row_names, col_names),
    index1   = FALSE,
    repr     = repr
  )
  
  attr(M, "flipped_snps") <- if (recode) col_names[flip] else character(0)
  return(M)
}