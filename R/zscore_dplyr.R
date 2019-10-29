
zscore_dplyr <- function(dat, sp_col, site_col, bioregion_col,
                   criterion = "significant", plot = FALSE){
  require(dplyr)
  require(ggplot2)

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

  if(!(criterion %in% c("significant", "top10"))){
    stop("criterion must be either 'significant' or 'top10'")
  }

  if(!is.logical(plot)){
    stop("plot must be a boolean.")
  }

  # Reassigning columns
  dat$sp <- dat[, sp_col]
  dat$site <- dat[, site_col]
  dat$bioregion <- dat[, bioregion_col]

  # Total number of sites
  n_site <- length(unique(dat$site))

  # Computing zscores => WRONG because absence are not taken into account
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
    as.data.frame()

  # Computing lambdas (connection between different bioregions)
  # Criterion: either zscore > 1.96 or top 10 zscores
  # Computing lambda
  if(criterion == "significant"){
    zscore_lambda <- zscore_sp %>%
      # NA if zscore inferior to 95% quantile of Gaussian distribution
      mutate(zscore = ifelse(zscore < 1.96, NA, zscore)) %>%
      # for each species: sum of significant zscore over the bioregions
      group_by(sp) %>%
      mutate(zscore_sum = sum(zscore, na.rm = TRUE),
             rho = zscore/zscore_sum) %>%
      ungroup() %>%
      as.data.frame()

  } else if(criterion == "top10"){
    # dat <- z_scores # for example
    zscore_lambda <- zscore_sp %>%
      # NA if zscore inferior to 90% quantile
      mutate(zscore = ifelse(zscore < as.numeric(quantile(zscore, 0.9,
                                                          na.rm = TRUE)),
                             NA, zscore)) %>%
      # for each species: sum of significant zscore over the bioregions
      group_by(sp) %>%
      mutate(zscore_sum = sum(zscore, na.rm = TRUE),
             rho = zscore/zscore_sum) %>%
      ungroup() %>%
      as.data.frame()
  }

  dat2 <- zscore_lambda[complete.cases(zscore_lambda), ]
  # Computing lambda indices
  lambda <- c()
  for(i in 1:length(unique(dat$bioregion))){
    focal_bioregion <- as.character(unique(dat$bioregion)[i])
    list_sp <- dat2[which(dat2$bioregion == unique(dat2$bioregion)[i]),
                    "sp"]
    list_sp <- as.character(list_sp)

    tmp <- dat %>%
      filter(sp %in% list_sp) %>%
      group_by(bioregion) %>%
      mutate(sum_rho = sum(rho, na.rm = TRUE)) %>%
      ungroup() %>%
      distinct(bioregion, .keep_all = TRUE) %>%
      mutate(focal_bioregion = focal_bioregion) %>%
      as.data.frame()

    lambda <- rbind(lambda,
                    tmp[, c("focal_bioregion", "bioregion", "sum_rho")])
  }

  # Remove duplicates per sp and site column from zscore
  zscore_sp <- zscore_sp %>%
    distinct(sp, .keep_all = TRUE) %>%
    dplyr::select(-site)

  if(plot == TRUE){
    res_plot <- ggplot(lambda, aes(focal_bioregion, sum_rho)) +
      geom_bar(aes(fill = as.factor(bioregion)), stat = "identity") +
      scale_fill_viridis_d("Bioregions") +
      labs(title = "Interaction between bioregions",
           x = "Bioregion", y = "Sum of contributions (%)") +
      theme_classic()
    return(list(zscore_sp, lambda, res_plot))
  } else{
    return(list(zscore_sp, lambda))
  }
}
