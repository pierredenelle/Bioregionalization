
contingency <- function(dat, site, sp, ab = NULL, binary = TRUE){
  require(DescTools)

  if(!is.data.frame(dat)){
    stop("dat must be a data.frame with columns sp and site.")
  }

  if(!is.character(site)){
    stop("site must be the column name of dat describing the sites")
  }

  if(!is.character(sp)){
    stop("sp must be the column name of dat describing the species.")
  }

  if(!is.null(ab) & !is.character(ab)){
    stop("ab must be the column name of dat describing the abundances
         of species")
  }

  if(is.null(ab) & binary == FALSE){
    warning("Without column abundances, contingency table will only get binary
        values.")
  }

  if(!is.logical(binary)){
    stop("binary must be a boolean.")
  }

  # Rename columns
  colnames(dat)[colnames(dat) == sp] <- "sp"
  colnames(dat)[colnames(dat) == site] <- "site"
  if(!is.null(ab)){
    colnames(dat)[colnames(dat) == ab] <- "ab"
  } else{ # without abundances, species are just present in the sites
    dat$ab <- 1
  }

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
