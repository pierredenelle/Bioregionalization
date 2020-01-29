
all_maps <- function(dat, form = "tidy", site, sp, ab = NULL, binary = TRUE,
                     similarity = "simpson", network_algo = "both",
                     saving_directory,
                     bipartite_algo = "greedy", weight = FALSE,
                     clustering = TRUE, ward_method = "ward.D2",
                     optim_method = "firstSEmax", nstart = 25, B = 50,
                     K.max = 20){
  ## Controls ----
  if(!is.character(form)){
    stop("form designs the format of the input data. It is either tidy (data
       frame with replicated rows per site) or a contingency table.")
  }

  if(!(form %in% c("tidy", "contingency"))){
    stop("form designs the format of the input data. It is either tidy (data
       frame with replicated rows per site) or a contingency table.")
  }

  if(!is.data.frame(dat)){
    stop("dat must be a data.frame with columns sp and site.")
  }

  if(!is.character(site)){
    stop("site must be the column name of dat describing the sites")
  }

  if(!is.character(sp)){
    stop("sp must be the column name of dat describing the species.")
  }

  if(!is.null(ab) & !is.character(ab)){
    stop("ab must be the column name of dat describing the abundances
         of species")
  }

  if(is.null(ab) & binary == FALSE){
    warning("Without column abundances, contingency table will only get binary
        values.")
  }

  if(!is.logical(binary)){
    stop("binary must be a boolean.")
  }

  if(!(similarity %in% c("simpson", "jaccard", "sorensen", "whittaker",
                         "bray"))){
    stop("Similarity metric chosen is not available.
     Please chose among the followings:
         simpson, jaccard, sorensen, whittaker or bray")
  }

  if(!is.character(network_algo)){
    stop("network_algo designs the type of algorithm among the following
         choices: projected, bipartite or both.")
  }

  if(!(network_algo %in% c("both", "projected", "bipartite"))){
    stop("network_algo designs the type of algorithm among the following
         choices: projected, bipartite or both.")
  }

  if(!(is.character(saving_directory))){
    stop("saving_directory must be a path where the OSLOM .tp file containing
         the bioregions identified will be saved.")
  }

  if(!(bipartite_algo %in% c("greedy", "girvan", "walktrap", "louvain", "LPAwb"))){
    stop("Provided algorithm to compute modularity is not available.
     Please chose among the followings:
         greedy, girvan, walktrap, louvain, or LPAwb.")
  }

  if(!is.logical(weight)){
    stop("weight must be a boolean.")
  }

  if(!is.logical(clustering)){
    stop("clustering must be a boolean indicating whether clustering techniques
         should be computed.")
  }

  if(!(method %in% c("ward.D", "ward.D2", "single", "complete", "average",
                     "mcquitty", "median", "centroid"))){
    stop("Hierarchical clustering method chosen is not available.
     Please chose among the followings:
         ward.D, ward.D2, single, complete, average,
         mcquitty, median or centroid.")
  }

  if(!(optim_method %in% c("globalmax", "firstmax", "Tibs2001SEmax",
                           "firstSEmax", "globalSEmax."))){
    stop("Chosen gap statistic to determine the optimal number of cluster is
    not available.
     Please chose among the followings:
         globalmax, firstmax, Tibs2001SEmax, firstSEmax or globalSEmax.")
  }

  if(!(abs(nstart - round(nstart)) < .Machine$double.eps^0.5)){
    stop("nstart must be an integer determining the number of random centroids
         to start k-means analysis.")
  }

  if(!(abs(B - round(B)) < .Machine$double.eps^0.5)){
    stop("B must be an integer determining the number of Monte Carlo bootstrap
         samples.")
  }

  if(!is.numeric(K.max)){
    stop("K.max must be a numeric determining the maximum number of clusters
         to consider.")
  }

  if(form == "contingency" & K.max > nrow(dat)){
    stop("K.max should not be superior to the number of rows of the
         contincenty matrix.")
  }

  ## Workflow ----
  # Create contingency table if needed
  if(form == "tidy"){
    dat <- contingency(dat, site = site, sp = sp, ab = ab, binary = binary)
  }

  # Project network
  dat_proj <- project_network(dat, similarity = similarity)
  dat_proj <- dat_proj[, c("id1", "id2", similarity)]

  # List of results
  list_res <- list()

  if(network_algo %in% c("projected", "both")){
    # Run OSLOM
    run_oslom(sp_proj, n_runs = 5, t_param = 0.1, cp_param = 0.5,
              saving_directory = saving_directory)
    res <- readRDS(paste0(saving_directory, "/tp.rds"))
    oslom_res <- oslom_output(res, dat)

    list_res["oslom"] <- oslom_res
  }

  if(network_algo %in% c("bipartite", "both")){
    # Bipartite methods
    bip <- algo_bipartite(dat, algo = bipartite_algo, weight = weight)

    list_res["bipartite"] <- bip
  }

  if(clustering == TRUE){
    # Clustering methods
    ward_res <- ward_cluster(
      dat, method = ward_method, optim_method = optim_method,
      nstart = nstart, B = B, K.max = K.max)

    list_res["ward"] <- ward_res
  }

  # Contribution of species
  tmp <- left_join(dat, list_res["oslom"], by = "site")
  scores <- contribute(dat = tmp, sp_col = "sp", site_col = "site",
                       bioregion_col = "bioregion")

  # Links between modules
  lambda <- interact(input_network = "projected",
                     dat = scores, plot = TRUE, output_format = "matrix")

  # Pixels all the time together
  # Gather all the bioregionalizations
  all_bioregions <- dat %>%
    select(site, x, y, env) %>%
    distinct(site, .keep_all = TRUE) %>%
    left_join(list_res["oslom"], by = "site") %>% # add OSLOM
    rename(oslom = bioregion) %>%
    left_join(list_res["ward"], by = "site") %>% # add Ward
    rename(ward = cluster)

  # Test of comparison function
  all100 <- comparison(all_bioregions, bio_col = c(5:6))

  return(list_res, scores, lambda, all100)
}
