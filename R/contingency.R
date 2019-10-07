
contingency <- function(dat, sp, site, ab, binary = FALSE){
  require(DescTools)
  # Rename columns
  colnames(dat)[colnames(dat) == sp] <- "sp"
  colnames(dat)[colnames(dat) == site] <- "site"
  colnames(dat)[colnames(dat) == ab] <- "ab"

  # Create contingency table
  mat <- xtabs(ab ~ site + sp, data = dat)
  mat <- as.matrix.xtabs(mat)

  # Check for empty rows and columns
  mat <- mat[rowSums(mat) > 0, ]
  mat <- mat[, colSums(mat) > 0]

  # Conversion as binary matrix if binary == TRUE
  if(binary == TRUE){
    mat[mat > 0] <- 1
  }

  return(mat)
}
