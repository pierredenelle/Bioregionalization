
oslom_output <- function(dat, contingency_mat){

  if(!is.matrix(contingency_mat)){
    stop("Contingency table should be a matrix with sites as rows and
         species as columns.")
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
    bioregion_oslom[which(bioregion_oslom$id %in% doublon$id), ]

    # Remove duplicates (first bioregion is assigned)
    bioregion_oslom <- bioregion_oslom[!duplicated(bioregion_oslom$id), ]
  }

  warning("WARNING: the order of the rownames of the contingency table
    has to be identical before running OSLOM and when running this function!!")

  return(bioregion_oslom[, c("site", "bioregion")])
}
