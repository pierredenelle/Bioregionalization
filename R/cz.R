
cz <- function(dat, sp_col, site_col, bip, ab = NULL){

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

  if(!is.data.frame(bip)){
    stop("bip must be a data.frame with three columns containing the nodes,
         their modules and their bipartite category.")
  }

  if(ncol(bip) != 3){
    stop("bip must be a data.frame with three columns containing the nodes,
         their modules and their bipartite category.")
  }

  if(sort(colnames(bip)) != c("cat", "mod", "node")){
    stop("Columns of bip must be named cat, mod and node.")
  }

  require(dplyr)

  # Construct table of links
  link_dat <- dat[, c("site", "sp")] %>%
    left_join(bip[which(bip$cat == "site"), c("node", "mod")],
              by = c("site" = "node")) %>%
    rename(mod_site = mod) %>%
    left_join(bip[which(bip$cat == "sp"), c("node", "mod")],
              by = c("sp" = "node")) %>%
    rename(mod_sp = mod)

  # Compute C
  C_site <- c()
  dat_site <- bip[which(bip$cat == "site"), ]
  for(i in 1:nrow(dat_site)){
    tmp2 <- table(link_dat[which(link_dat$site == dat_site[i, "node"]),
                           "mod_sp"])
    C_site <- rbind(C_site,
                    c(as.character(dat_site[i, "node"]),
                      1 - sum((tmp2/sum(tmp2))^2), "site"))
  }
  C_sp <- c()
  dat_sp <- bip[which(bip$cat == "sp"), ]
  for(i in 1:nrow(dat_sp)){
    tmp2 <- table(link_dat[which(link_dat$sp == dat_sp[i, "node"]), "mod_site"])
    C_sp <- rbind(C_sp,
                  c(as.character(dat_sp[i, "node"]),
                    1 - sum((tmp2/sum(tmp2))^2), "sp"))
  }
  C_dat <- rbind(C_site, C_sp)
  C_dat <- as.data.frame(C_dat)
  colnames(C_dat) <- c("node", "C", "cat")

  bip <- left_join(bip, C_dat, by = c("node",  "cat"))
  bip$C <- as.numeric(as.character(bip$C))

  # Compute z
  bip$n_link_mod <- NA
  for(i in 1:nrow(bip)){
    if(bip[i, "cat"] == "site"){
      tmp2 <- link_dat[which(link_dat$site == bip[i, "node"]), ]
      bip[i, "n_link_mod"] <- nrow(tmp2[which(tmp2$mod_site == tmp2$mod_sp), ])
    } else{
      tmp2 <- link_dat[which(link_dat$sp == bip[i, "node"]), ]
      bip[i, "n_link_mod"] <- nrow(tmp2[which(tmp2$mod_site == tmp2$mod_sp), ])
    }
  }

  bip$mod <- as.character(bip$mod)

  mean_link_mod <- tapply(bip$n_link_mod, bip$mod, mean)
  mean_link_mod <- data.frame("mod" = names(mean_link_mod),
                              "mean_link_mod" = as.numeric(mean_link_mod))
  bip <- left_join(bip, mean_link_mod, by = "mod")

  sd_link_mod <- tapply(bip$n_link_mod, bip$mod, sd)
  sd_link_mod <- data.frame("mod" = names(sd_link_mod),
                            "sd_link_mod" = as.numeric(sd_link_mod))
  bip <- left_join(bip, sd_link_mod, by = "mod")

  bip$z <- (bip$n_link_mod - bip$mean_link_mod) / bip$sd_link_mod

  # Z with standardisation only for species
  dat_sp <- bip[which(bip$cat == "sp"), ]
  mean_link_mod_sp <- tapply(dat_sp$n_link_mod, dat_sp$mod, mean)
  mean_link_mod_sp <- data.frame("mod" = names(mean_link_mod_sp),
                                 "mean_link_mod_sp" = as.numeric(mean_link_mod_sp))
  dat_sp <- left_join(dat_sp, mean_link_mod_sp, by = "mod")

  sd_link_mod_sp <- tapply(dat_sp$n_link_mod, dat_sp$mod, sd)
  sd_link_mod_sp <- data.frame("mod" = names(sd_link_mod_sp),
                               "sd_link_mod_sp" = as.numeric(sd_link_mod_sp))
  dat_sp <- left_join(dat_sp, sd_link_mod_sp, by = "mod")
  dat_sp$z <- (dat_sp$n_link_mod - dat_sp$mean_link_mod_sp) / dat_sp$sd_link_mod_sp

  return(list(bip, dat_sp))

}
