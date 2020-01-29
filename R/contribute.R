
contribute <- function(dat, sp_col, site_col, bioregion_col,
                       bioregion_sp_col = NULL, ab = NULL){

  # Controls and initial steps ------------------------------------------------
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

  # Reassigning columns
  dat$sp <- dat[, sp_col]
  dat$site <- dat[, site_col]
  dat$bio_site <- dat[, bioregion_col]

  if(!is.null(bioregion_sp_col)){
    if(!is.character(bioregion_col)){
      stop("bioregion_sp_col must be a character string corresponding to the
      column with the bioregions of species.")
    }
  } else{
    warning("No bioregion provided for species. Each species will be assigned
            to the bioregion where it occurs the most.")
    # Assign bioregion to species
    sp_bio <- dat %>%
      group_by(sp, bioregion) %>%
      summarise(c = n()) %>%
      filter(row_number(desc(c)) == 1) %>%
      rename(bio_sp = bioregion) %>%
      select(sp, bio_sp)
    # Merge with site-species data frame
    dat <- left_join(dat, sp_bio, by = "sp")
  }

  require(dplyr)

  # C -------------------------------------------------------------------------
  # Guimera et Amaral 2005 Nature
  C_site <- c()
  dat_site <- dat[!duplicated(dat$site), ]
  for(i in 1:nrow(dat_site)){
    tmp <- table(dat[which(dat$site == dat_site[i, "site"]), "bio_sp"])
    C_site <- rbind(C_site,
                    c(as.character(dat_site[i, "site"]),
                      1 - sum((tmp/sum(tmp))^2)))
  }
  C_site <- as.data.frame(C_site)
  colnames(C_site) <- c("site", "C_site")

  C_sp <- c()
  dat_sp <- dat[!duplicated(dat$sp), ]
  for(i in 1:nrow(dat_sp)){
    tmp <- table(dat[which(dat$sp == dat_sp[i, "sp"]), "bio_site"])
    C_sp <- rbind(C_sp,
                  c(as.character(dat_sp[i, "sp"]),
                    1 - sum((tmp/sum(tmp))^2)))
  }
  C_sp <- as.data.frame(C_sp)
  colnames(C_sp) <- c("sp", "C_sp")

  # Merge with site-species data frame
  dat <- left_join(dat, C_site, by = "site")
  dat <- left_join(dat, C_sp, by = "sp")

  # z -------------------------------------------------------------------------
  dat_site$n_link_mod <- NA
  for(i in 1:nrow(dat_site)){
    tmp <- dat[which(dat$site == dat_site[i, "site"]), ]
    dat_site[i, "n_link_mod"] <- nrow(tmp[which(tmp$bio_site == tmp$bio_sp), ])
  }

  dat_sp$n_link_mod <- NA
  for(i in 1:nrow(dat_sp)){
    tmp <- dat[which(dat$sp == dat_sp[i, "sp"]), ]
    dat_sp[i, "n_link_mod"] <- nrow(tmp[which(tmp$bio_sp == tmp$bio_sp), ])
  }

  dat_site$bio_site <- as.character(dat_site$bio_site)
  dat_sp$bio_sp <- as.character(dat_sp$bio_sp)

  dat_site <- dat_site %>%
    group_by(bio_site) %>%
    mutate(mean_link_mod = mean(n_link_mod),
           sd_link_mod = sd(n_link_mod),
           z_site = (n_link_mod - mean_link_mod)/sd_link_mod)

  dat_sp <- dat_sp %>%
    group_by(bio_sp) %>%
    mutate(mean_link_mod = mean(n_link_mod),
           sd_link_mod = sd(n_link_mod),
           z_sp = (n_link_mod - mean_link_mod)/sd_link_mod)

  # Merge with site-species data frame
  dat <- left_join(dat, dat_site[, c("site", "bio_site", "z_site")],
                   by = c("site", "bio_site"))
  dat <- left_join(dat, dat_sp[, c("sp", "bio_sp", "z_sp")],
                   by = c("sp", "bio_sp"))

  # rho -----------------------------------------------------------------------
  # See Lenormand et al. 2019 Ecology and Evolution
  # Reassigning columns
  dat$sp <- dat[, sp_col]
  dat$site <- dat[, site_col]
  dat$bioregion <- dat[, bioregion_col]

  # Number of cells per bioregions
  clust <- table(dat[!duplicated(dat[, "site"]), "bio_site"])
  clust <- data.frame(bio_site = names(clust), ncell = as.numeric(clust))

  n_bioregion <- length(unique(dat[complete.cases(dat), "bio_site"]))
  n_site <- length(unique(dat[complete.cases(dat), "site"]))

  # Compute nij, ni & nj
  agg <- aggregate(dat$site, list(dat$sp, dat$bio_site), length)
  colnames(agg) <- c("sp", "bio_site", "nij")

  # Number of cells in bioregion j where the specie i is present
  nij <- xtabs(nij ~ sp + bio_site, data = agg)
  # Number of cells where the specie i is present
  ni <- replicate(ncol(nij), apply(nij, 1, sum))
  # Number of cells in the bioregion j
  nj <- t(replicate(nrow(nij), clust[, "ncell"]))

  # Rhoij
  num <- nij-((ni*nj)/n_site)
  den <- sqrt((n_site-nj)/(n_site-1)*(1-(nj/n_site))*((ni*nj)/n_site))
  rhoij <- num/den

  # Conversion to data frame
  rhoij <- as.data.frame.matrix(rhoij) # species in rows & bioregions in columns
  # data.frame format for rhoij
  rhoij <- as.data.frame(as.table(as.matrix(rhoij)))
  colnames(rhoij) <- c("sp", "bio_site", "rho")

  # Merge with site-species data frame
  dat <- left_join(dat, rhoij, by = c("sp", "bio_site"))

  # Return data frame for species ---------------------------------------------
  return(dat)

}
