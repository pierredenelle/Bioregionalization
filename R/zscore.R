
zscore <- function(dat, sp_col, site_col, bioregion_col){
  require(dplyr)

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

  # Reassigning columns
  dat$sp <- dat[, sp_col]
  dat$site <- dat[, site_col]
  dat$bioregion <- dat[, bioregion_col]

  # Total number of sites
  n_site <- length(unique(dat$site))

  # Computing zscores
  zscore_sp <- dat %>%
    # occurrences of species
    group_by(sp) %>%
    mutate(n_i = n()) %>%
    ungroup() %>%
    # number of sites of each bioregion
    group_by(bioregion) %>%
    mutate(n_j = length(unique(site))) %>%
    ungroup() %>%
    # zscore
    group_by(bioregion, sp) %>%
    mutate(n_ij = n(), # occurrence of sp i in bioregion j
           zscore = (n_ij - n_i*n_j/n_site)/
             sqrt((n_site - n_j)/(n_site-1)*(1-n_j/n_site)*n_i*n_j/n_site)) %>%
    distinct(sp, .keep_all = TRUE) %>% # removal of duplicates per sp
    dplyr::select(-site) %>% # removal of pixel column
    as.data.frame()

  return(zscore_sp)
}
