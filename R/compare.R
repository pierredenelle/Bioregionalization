
comparison <- function(dat){

  if(!is.data.frame(dat)){
    stop("Input must be a data.frame with a column containing the sites and
         columns containing different partitions.")
  }

  if(!is.character(sp_col)){
    stop("sp_col must be a character string corresponding to the column
    with species.")
  }

  # Identifying common pixels through different methods
  group1 <- all_bioregions[, "site"]

  list1_oslom <-
    all_bioregions[which(all_bioregions$oslom == unique(all_bioregions$oslom)[1]),
                   "site"]
  list1_ward <- all_bioregions[which(all_bioregions$ward == unique(all_bioregions$ward)[1]),
                               "site"]
  list1_greedy <- all_bioregions[which(all_bioregions$greedy == unique(all_bioregions$greedy)[1]),
                                 "site"]
  list1_lpawb <- all_bioregions[which(all_bioregions$lpawb == unique(all_bioregions$lpawb)[1]),
                                "site"]
  list1_infomap <- all_bioregions[which(all_bioregions$infomap == unique(all_bioregions$infomap)[1]),
                                  "site"]

  list1_combined <- list1_oslom[list1_oslom %in% list1_ward]
  list1_combined <- list1_combined[list1_combined %in% list1_greedy]
  list1_combined <- list1_combined[list1_combined %in% list1_lpawb]
  list1_combined <- list1_combined[list1_combined %in% list1_infomap]

  all_bioregions$combined <- NA
  all_bioregions[which(all_bioregions$site %in% list1_combined), "combined"] <- "1"

  # Output: data frame with pixels grouped together through all the methods
    return(combined_groups)
}
