
lambda <- function(dat, sp_col, zscore_col, bioregion_col,
                   criterion = "significant", plot = FALSE){
  require(dplyr)
  require(ggplot2)

  if(!is.character(sp_col)){
    stop("sp_col must be a character string corresponding to the column
    with species.")
  }

  if(!is.character(zscore_col)){
    stop("zscore_col must be a character string corresponding to the column
    with zscores per species/biroegion combination.")
  }

  if(!is.character(bioregion_col)){
    stop("bioregion_col must be a character string corresponding to the column
    with bioregions.")
  }

  if(!is.numeric(dat[, zscore_col])){
    stop("zscore column must contain numeric values corresponding to
         the zscores of species (see zscore function).")
  }

  if(!(criterion %in% c("significant", "top10"))){
    stop("criterion must be either 'significant' or 'top10'")
  }

  if(!is.logical(plot)){
    stop("plot must be a boolean.")
  }

  # Reassigning columns
  dat$sp <- dat[, sp_col]
  dat$zscore_col <- dat[, zscore_col]
  dat$bioregion <- dat[, bioregion_col]

  # Criterion: either zscore>1.96 or top 10 zscores
  # Computing lambda
  if(criterion == "significant"){
    dat <- dat %>%
      # NA if zscore inferior to 90% quantile
      mutate(zscore = ifelse(zscore < 1.96, NA, zscore)) %>%
      group_by(bioregion) %>%
      mutate(zscore_sum = sum(zscore, na.rm = TRUE),
             rho = zscore/zscore_sum) %>%
      ungroup() %>%
      as.data.frame()

  } else if(criterion == "top10"){
    # dat <- z_scores # for example
    dat <- dat %>%
      # NA if zscore inferior to 90% quantile
      mutate(zscore = ifelse(zscore < as.numeric(quantile(zscore, 0.9,
                                                          na.rm = TRUE)),
                             NA, zscore)) %>%
      group_by(bioregion) %>%
      mutate(zscore_sum = sum(zscore, na.rm = TRUE),
             rho = zscore/zscore_sum) %>%
      ungroup() %>%
      as.data.frame()
  }

  dat2 <- dat[complete.cases(dat), ]
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

  if(plot == TRUE){
    res_plot <- ggplot(lambda, aes(focal_bioregion, sum_rho)) +
      geom_bar(aes(fill = as.factor(bioregion)), stat = "identity") +
      scale_fill_viridis_d("Bioregions") +
      labs(title = "Interaction between bioregions",
           x = "Bioregion", y = "Sum of contributions (%)") +
      theme_classic()
    return(list(lambda, res_plot))
  } else{
    return(lambda)
  }
}
