
oslom_output <- function(dat, contingency_mat, doublon_removal = "first"){

  require(dplyr)

  if(!is.matrix(contingency_mat)){
    stop("Contingency table should be a matrix with sites as rows and
         species as columns.")
  }

  if(!(doublon_removal %in% c("first", ""))){
    stop("doublon_removal must be a character string equal to 'first' or
         'max'.")
  }

  if(is.null(rownames(contingency_mat))){
    stop("contingency_mat should have the sites set as rownames.")
  }

  # Convert dat into list
  bioregion_list <- list()
  length(bioregion_list) <- length(dat)/2
  for(k in 1:length(dat)){
    if((k/2-trunc(k/2)) == 0){ # loop over module pixels
      bioregion_list[[(k/2)]] <-
        as.numeric(as.matrix(strsplit(dat[k], split = " ")[[1]]))
    }
  }

  # Convert bioregion_list into data.frame
  names(bioregion_list) <- as.character(seq_along(bioregion_list))
  bioregion_oslom <- data.frame(
    id_oslom = unlist(bioregion_list),
    bioregion = rep(names(bioregion_list), sapply(bioregion_list, length)))

  # Get original site from the relevant contingency table
  # WARNING: the order of the rownames of the contingency table
  # has to be identical before running OSLOM and when running this function
  bioregion_oslom$site <- rownames(contingency_mat)[bioregion_oslom$id_oslom]

  # If some doublons are present, attribute first bioregion
  if(length(bioregion_list) > 1 &
     length(unique(table(bioregion_oslom$id))) > 1){
    # Identify doublons and get the first one
    doublon <- table(bioregion_oslom$id)[table(bioregion_oslom$id) > 1]
    doublon <- data.frame(id = names(doublon),
                          nb_bioregion = as.numeric(doublon))

    if(doublon_removal == "first"){
      # Remove duplicates (first bioregion is assigned)
      bioregion_oslom <- bioregion_oslom[!duplicated(bioregion_oslom$id), ]
    } else if(doublon_removal == "max"){
      # Maximum similarity => get score for each species
      # then assign bioregion for which the score is maximal
      # bioregion_oslom[which(bioregion_oslom$id %in% doublon$id), ]
      score_sp <- c()
      sp_mat_without_doublon <-
        sp_mat[-as.numeric(as.character(doublon$id)), ]
      for(i in 1:ncol(sp_mat)){
        pres_com <-
          names(sp_mat_without_doublon[, i][sp_mat_without_doublon[, i] > 0])
        # Bioregion of these sites
        max_bioregion <- table(
          bioregion_oslom[which(bioregion_oslom$site %in% pres_com),
                          "bioregion"])
        max_bioregion <- names(which.max(max_bioregion))
        score_sp <- rbind(score_sp,
                          data.frame(sp = colnames(sp_mat)[i],
                                     max_bioregion = max_bioregion))
      }

      doublon$bioregion <- NA
      for(i in 1:length(doublon$id)){
        i_sites <- sp_mat[as.numeric(as.character(doublon[i, "id"])), ]
        i_sites <- i_sites[i_sites > 0]
        choice <- score_sp[which(score_sp$sp %in% names(i_sites)),
                           "max_bioregion"]
        choice <- table(choice)
        choice <- names(which.max(choice))
        doublon[i, "bioregion"] <- choice
      }
      bioregion_oslom$id_oslom <- as.character(bioregion_oslom$id_oslom)
      colnames(doublon)[colnames(doublon) == "id"] <- "id_oslom"
      doublon <- left_join(doublon, bioregion_oslom[, c("id_oslom", "site")],
                           by = c("id_oslom"))

      bioregion_oslom <-
        rbind(
          bioregion_oslom[which(bioregion_oslom$id %in% doublon$id), ],
          doublon[, c("id_oslom", "bioregion", "site")])
    }
    # CREATE FUNCTION REASSIGN => also Ward
  }

  warning("WARNING: the order of the rownames of the contingency table
    has to be identical before running OSLOM and when running this function!!")

  return(bioregion_oslom[, c("site", "bioregion")])
}
