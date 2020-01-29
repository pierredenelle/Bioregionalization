
interact <- function(input_network, dat, plot = FALSE,
                     output_format = "matrix"){

  if(!is.character(input_network)){
    stop("input_network must be a character string equal to 'projected' or
         'bipartite'. It determines the format of the input.")
  }

  if(input_network == "projected"){
    if(!is.matrix(dat)){
      stop("When input_network is 'projected', input must be a data.frame
      with each row indicating the contribution of a species in a given
           bioregion.")
    }

    if(is.null(rownames(dat))){
      stop("When input_network is 'projected', input matrix must have species
           names as rownames.")
    }

    if(is.null(colnames(dat))){
      stop("When input_network is 'projected', input matrix must have bioregion
           names as colnames.")
    }
  } else if(input_network == "bipartite"){
    if(!is.data.frame(dat)){
      stop("When input_network is 'bipartite',
      dat must be a data.frame with four columns containing the links between
           the two categories of nodes and their respective modules.")
    }

    if(ncol(dat) != 4){
      stop("When input_network is 'bipartite',
      dat must be a data.frame with four columns containing the links between
           the two categories of nodes and their respective modules.")
    }

    if(sort(colnames(dat)) != c("mod_site", "mod_sp", "site", "sp")){
      stop("When input_network is 'bipartite', columns of dat must be named
           mod_site, mod_sp, site and sp.")
    }
  } else{
    stop("input_network must be a character string equal to 'projected' or
         'bipartite'. It determines the format of the input.")
  }

  if(!is.logical(plot)){
    stop("plot argument should be a boolean determining whether the barplot of
         interactions between bioregions should be saved.")
  }

  if(!is.character(output_format)){
    stop("output_format must be a character string equal to 'matrix' or
         'data.frame'. It determines the format of the output.")
  }

  if(!(output_format %in% c("matrix", "dataframe"))){
    stop("output_format must be a character string equal to 'matrix' or
         'data.frame'. It determines the format of the output.")
  }

  if(input_network == "projected"){
    # Lambda
    rhoijp <- dat
    # NA if zscore inferior to 95% quantile of Gaussian distribution
    # rhoijp[rhoijp < 1.96] <- NA
    # NA if zscore inferior to mean of the contribution per bioregion
    rhoijp[sweep(rhoijp, 2, colMeans(rhoijp), "<")] <- NA
    # NA if zscore negative
    rhoijp[rhoijp < 0] <- NA

    # for each species: sum of significant rhos over the bioregions
    # dim(dat) = number of species (rows) and number of bioregions (columns)
    rhoijp <- rhoijp/rowSums(rhoijp, na.rm = TRUE)

    lambda <- NULL
    # Bioregions in columns of input matrix
    n_bioregion <- ncol(dat)

    for(k in 1:n_bioregion){ # loop over the bioregions
      rhoijpk <- rhoijp[!is.na(rhoijp[, k]), ]
      # Control for cases where only one species if assigned to one module
      if(!is.null(dim(rhoijpk))){
        lambda <- rbind(
          lambda,
          100*apply(rhoijpk, 2, sum, na.rm = TRUE)/nrow(rhoijpk))
      } else{
        rhoijpk <- t(as.matrix(rhoijpk))
        lambda <- rbind(
          lambda,
          100*apply(rhoijpk, 2, sum, na.rm = TRUE)/nrow(rhoijpk))
      }
    }

    # Save results
    rownames(lambda) <- colnames(dat)
    colnames(lambda) <- rownames(lambda)

    # Data frame format for plot
    lambda_plot <- as.data.frame(as.table(lambda))
    colnames(lambda_plot) <- c("bioregion", "link_bioregion", "lambda")

    if(output_format == "dataframe"){
      lambda <- lambda_plot
      colnames(lambda) <- c("focal_bioregion", "bioregion", "lambda")
    }

    # Output: the matrix of bioregion relationships lambda
    if(plot == TRUE){
      # Barplot
      res_plot <- ggplot(lambda_plot, aes(bioregion, lambda)) +
        geom_bar(aes(fill = as.factor(link_bioregion)), stat = "identity") +
        scale_fill_viridis_d("Bioregions") +
        labs(title = "Interaction between bioregions",
             x = "Bioregion", y = "Sum of contributions (%)") +
        theme_classic() +
        theme(panel.border = element_rect(fill = NA))

      return(list(lambda = lambda, res_plot))
    } else{
      return(lambda = lambda)
    }

  } else if(input_network == "bipartite"){
    # Percentage of species from the different modules per module
    lambda <- table(dat$mod_site, dat$mod_sp)
    lambda <- 100*lambda/rowSums(lambda)

    # Data frame format for plot
    lambda_plot <- as.data.frame(as.table(lambda))
    colnames(lambda_plot) <- c("mod_site", "mod_sp", "perc")

    if(output_format == "dataframe"){
      lambda <- lambda_plot
    }

    # Output: the matrix of bioregion relationships lambda
    if(plot == TRUE){
      # Barplot
      res_plot <- ggplot(lambda_plot, aes(mod_site, perc)) +
        geom_bar(aes(fill = as.factor(mod_sp)), stat = "identity") +
        scale_fill_viridis_d("Bioregions") +
        labs(title = "Interaction between bioregions",
             x = "Bioregion", y = "Sum of contributions (%)") +
        theme_classic() +
        theme(panel.border = element_rect(fill = NA))

      return(list(lambda = lambda, res_plot))
    } else{
      return(lambda = lambda)
    }
  }
}
