
zscore <- function(dat, sp_col, site_col, bioregion_col,
                   output_format = "matrix"){

  if(!is.data.frame(dat)){
    stop("Input must be a data.frame with each row indicating the presence
         of a species in a given site. The associated bioregion must be present
         for each site.")
  }

  if(!is.character(sp_col)){
    stop("sp_col must be a character string corresponding to the column
    with species.")
  }

  if(!is.character(site_col)){
    stop("site_col must be a character string corresponding to the column
    with sites")
  }

  if(!is.character(bioregion_col)){
    stop("bioregion_col must be a character string corresponding to the column
    with bioregions.")
  }

  if(!is.character(output_format)){
    stop("output_format must be a character string equal to 'matrix' or
         'data.frame'. It determines the format of the output.")
  }

  if(!(output_format %in% c("matrix", "dataframe"))){
    stop("output_format must be a character string equal to 'matrix' or
         'data.frame'. It determines the format of the output.")
  }

  # Reassigning columns
  dat$sp <- dat[, sp_col]
  dat$site <- dat[, site_col]
  dat$bioregion <- dat[, bioregion_col]

  # Number of cells per bioregions
  clust <- table(dat[!duplicated(dat[, "site"]), "bioregion"])
  clust <- data.frame(bioregion = names(clust), ncell = as.numeric(clust))

  n_bioregion <- length(unique(dat[complete.cases(dat), "bioregion"]))
  n_site <- length(unique(dat[complete.cases(dat), "site"]))

  # Compute nij, ni & nj
  agg <- aggregate(dat$site, list(dat$sp, dat$bioregion), length)
  colnames(agg) <- c("sp", "bioregion", "nij")

  # Number of cells in bioregion j where the specie i is present
  nij <- xtabs(nij ~ sp + bioregion, data = agg)
  # Number of cells where the specie i is present
  ni <- replicate(ncol(nij), apply(nij, 1, sum))
  # Number of cells in the bioregion j
  nj <- t(replicate(nrow(nij), clust[, "ncell"]))

  # Rhoij
  num <- nij-((ni*nj)/n_site)
  den <- sqrt((n_site-nj)/(n_site-1)*(1-(nj/n_site))*((ni*nj)/n_site))
  rhoij <- num/den

  if(output_format == "dataframe"){
    rhoij <- as.data.frame.matrix(rhoij) # species in rows & bioregions in columns
    # data.frame format for rhoij
    rhoij <- as.data.frame(as.table(as.matrix(rhoij)))
    colnames(rhoij) <- c("sp", "bioregion", "zscore")
  }

  # Output: the test-value matrix rhoij
  return(zscore_sp = rhoij)
}
