
comparison <- function(dat, bio_col, thres = 10){

  if(!is.data.frame(dat)){
    stop("Input must be a data.frame with a column containing the sites and
         columns containing different partitions.")
  }

  if(!("site" %in% colnames(dat))){
    stop("dat must contains a column 'site' with all the sites.")
  }

  if(!is.numeric(bio_col) & bio_col > ncol(dat)){
    stop("bio_col must be a numeric indicating which columns of dat contains
         the different bioregions.")
  }

  if(!is.numeric(thres)){
    stop("thres must be a numeric determining the minimum number of pixels
         per pair.")
  }

  # Get matrix with number of bioregions in common per pair of pixels
  list_mat_pair <- list()
  for(i in 1:length(bio_col)){
    list_mat_pair[[i]] <-
      matrix(outer(dat[, bio_col[i]], dat[, bio_col[i]], "=="),
             nrow = nrow(dat), dimnames = list(dat$site, dat$site))
  }
  # Sum all matrices
  list_mat_pair <- Reduce('+', list_mat_pair)

  # Get percentage
  list_mat_pair <- 100*list_mat_pair/length(bio_col)

  # Conversion to data frame
  list_df_pair <- as.data.frame(as.table(list_mat_pair))
  colnames(list_df_pair) <- c("id1", "id2", "perc")

  # Remove false pairs of sites
  list_df_pair <- list_df_pair[which(list_df_pair$id1 != list_df_pair$id2), ]

  # Extract pixels together 100% of the time
  all100 <- list_df_pair[which(list_df_pair$perc == 100), ]

  # Add bioregion (whatever algorithm since all pixels are grouped together)
  all100 <- left_join(all100,
                      all_bioregions[, c("site",
                                         colnames(all_bioregions)[bio_col])],
                      by = c("id1" = "site"))
  colnames(all100) <- c("id1", "id2", "perc",
                        paste0(colnames(all_bioregions)[bio_col], "_id1"))

  all100 <- left_join(all100,
                      all_bioregions[, c("site",
                                         colnames(all_bioregions)[bio_col])],
                      by = c("id2" = "site"))
  colnames(all100) <- c("id1", "id2", "perc",
                        paste0(colnames(all_bioregions)[bio_col], "_id1"),
                        paste0(colnames(all_bioregions)[bio_col], "_id2"))

  bioregion_cols <- c(paste0(colnames(all_bioregions)[bio_col], "_id1"),
                      paste0(colnames(all_bioregions)[bio_col], "_id2"))

  # all100$pair_bio <- paste(all100$bio_id1, all100$bio_id2, sep = "_")

  # all100$pair_bio <- tidyr::unite_(data,
  #                                  paste(bioregion_cols, collapse = "_"),
  #                                  bioregion_cols)
  all100$pair_bio <- apply(all100[, bioregion_cols] , 1, paste, collapse = "_")

  # Remove pairs of pixels under a threshold of 10 pixels
  thresh <- names(table(all100$pair_bio)[table(all100$pair_bio) > thres])

  all100 <- all100[which(all100$pair_bio %in% thresh), ]

  # Output: data frame with pixels grouped together through all the methods
  return(all100)
}
