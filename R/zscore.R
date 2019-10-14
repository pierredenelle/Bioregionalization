
zscore <- function(dat, sp_col, site_col, bioregion_col){
  require(dplyr)
  dat$sp <- dat[, sp_col]
  dat$site <- dat[, site_col]
  dat$bioregion <- dat[, bioregion_col]

  n_site <- length(unique(dat$site))

  zscore_sp <- dat %>%
    group_by(bioregion, sp) %>%
    mutate(n_i = n()) %>% # occurrence of species
    ungroup() %>%
    group_by(bioregion) %>%
    mutate(n_j = length(unique(site))) %>% # number of sites of each bioregion
    ungroup() %>%
    group_by(bioregion, sp) %>%
    mutate(n_ij = n(), # occurrence of sp i in bioregion j
           zscore = (n_ij - n_i*n_j/n_site)/
             sqrt((n_site - n_j)/(n_site-1)*(1-n_j/n_site)*n_i*n_j/n_site)) %>%
    distinct(sp, .keep_all = TRUE) %>% # removal of duplicates per sp
    dplyr::select(-site) %>% # removal of pixel column
    as.data.frame()

  return(zscore_sp)
}
