
##-----------------------------------------------------------------------------
# Function to compute Cz coefficients (Guimera et Amaral 2005)
# Needs two data frames as inputs:

# link_dat is a data.frame of the links between nodes containing four columns:
# site, sp, mod_sp, mod_site which stand for the first node category, the
# second node category, module of the first node and module of the second one

# dat contains three columns: node, mod and cat which respectively stand for
# the name of the node, its module and its bipartite category

# column of abundance can be added to weight the calculation
# (not implemented yet)

# The function returns a list containing two data.frame which correspond to
# dat data.frame with extra-columns Cz, z being only standardized per species
# in the second data.frame of the list

##-----------------------------------------------------------------------------

cz <- function(link_dat, dat, ab = NULL){

  if(!is.data.frame(link_dat)){
    stop("link_dat must be a data.frame with four columns containing the links
    between the two categories of nodes and their respective modules.")
  }

  if(ncol(link_dat) != 4){
    stop("link_dat must be a data.frame with four columns containing the links
    between the two categories of nodes and their respective modules.")
  }

  if(sort(colnames(link_dat)) != c("mod_site", "mod_sp", "site", "sp")){
    stop("Columns of dat must be named mod_site, mod_sp, site and sp.")
  }

  if(!is.data.frame(dat)){
    stop("dat must be a data.frame with three columns containing the nodes,
         their modules and their bipartite category.")
  }

  if(ncol(dat) != 3){
    stop("dat must be a data.frame with three columns containing the nodes,
         their modules and their bipartite category.")
  }

  if(sort(colnames(dat)) != c("cat", "mod", "node")){
    stop("Columns of dat must be named cat, mod and node.")
  }

  require(dplyr)

  # Compute C
  C_site <- c()
  dat_site <- dat[which(dat$cat == "site"), ]
  for(i in 1:nrow(dat_site)){
    tmp2 <- table(link_dat[which(link_dat$site == dat_site[i, "node"]),
                           "mod_sp"])
    C_site <- rbind(C_site,
                    c(as.character(dat_site[i, "node"]),
                      1 - sum((tmp2/sum(tmp2))^2), "site"))
  }
  C_sp <- c()
  dat_sp <- dat[which(dat$cat == "sp"), ]
  for(i in 1:nrow(dat_sp)){
    tmp2 <- table(link_dat[which(link_dat$sp == dat_sp[i, "node"]), "mod_site"])
    C_sp <- rbind(C_sp,
                  c(as.character(dat_sp[i, "node"]),
                    1 - sum((tmp2/sum(tmp2))^2), "sp"))
  }
  C_dat <- rbind(C_site, C_sp)
  C_dat <- as.data.frame(C_dat)
  colnames(C_dat) <- c("node", "C", "cat")

  dat <- left_join(dat, C_dat, by = c("node",  "cat"))
  dat$C <- as.numeric(as.character(dat$C))

  # Compute z
  dat$n_link_mod <- NA
  for(i in 1:nrow(dat)){
    if(dat[i, "cat"] == "site"){
      tmp2 <- link_dat[which(link_dat$site == dat[i, "node"]), ]
      dat[i, "n_link_mod"] <- nrow(tmp2[which(tmp2$mod_site == tmp2$mod_sp), ])
    } else{
      tmp2 <- link_dat[which(link_dat$sp == dat[i, "node"]), ]
      dat[i, "n_link_mod"] <- nrow(tmp2[which(tmp2$mod_site == tmp2$mod_sp), ])
    }
  }

  dat$mod <- as.character(dat$mod)

  mean_link_mod <- tapply(dat$n_link_mod, dat$mod, mean)
  mean_link_mod <- data.frame("mod" = names(mean_link_mod),
                              "mean_link_mod" = as.numeric(mean_link_mod))
  dat <- left_join(dat, mean_link_mod, by = "mod")

  sd_link_mod <- tapply(dat$n_link_mod, dat$mod, sd)
  sd_link_mod <- data.frame("mod" = names(sd_link_mod),
                            "sd_link_mod" = as.numeric(sd_link_mod))
  dat <- left_join(dat, sd_link_mod, by = "mod")

  dat$z <- (dat$n_link_mod - dat$mean_link_mod) / dat$sd_link_mod

  # Z with standardisation only for species
  dat_sp <- dat[which(dat$cat == "sp"), ]
  mean_link_mod_sp <- tapply(dat_sp$n_link_mod, dat_sp$mod, mean)
  mean_link_mod_sp <- data.frame("mod" = names(mean_link_mod_sp),
                                 "mean_link_mod_sp" = as.numeric(mean_link_mod_sp))
  dat_sp <- left_join(dat_sp, mean_link_mod_sp, by = "mod")

  sd_link_mod_sp <- tapply(dat_sp$n_link_mod, dat_sp$mod, sd)
  sd_link_mod_sp <- data.frame("mod" = names(sd_link_mod_sp),
                               "sd_link_mod_sp" = as.numeric(sd_link_mod_sp))
  dat_sp <- left_join(dat_sp, sd_link_mod_sp, by = "mod")
  dat_sp$z <- (dat_sp$n_link_mod - dat_sp$mean_link_mod_sp) / dat_sp$sd_link_mod_sp

  return(list(dat, dat_sp))

}
