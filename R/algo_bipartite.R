
algo_bipartite <- function(dat, algo = "greedy", weight = FALSE){

  if(!is.matrix(dat)){
    stop("Input dat should be a matrix with sites as rows and species
         as columns.")
  }

  if(!(algo %in% c("greedy", "girvan", "walktrap", "louvain", "LPAwb", "infomap"))){
    stop("Provided algorithm to compute modularity is not available.
     Please chose among the followings:
         greedy, girvan, walktrap, louvain, or LPAwb.")
  }

  if(!is.logical(weight)){
    stop("weight must be a boolean.")
  }

  require(igraph)
  # source("scripts/Beckett_LPAwb_function.R") => store in R function

  if(algo == "LPAwb"){ # Beckett modularity
    dat <- as.matrix(dat)
    # Find labels and weighted modularity using LPAwb+
    network_mod <- Bioregionalization::LPA_wb_plus(dat)

    # Conversion into data.frame with node, category and module
    network_lab <- data.frame(node = c(rownames(dat), colnames(dat)),
                              module = c(as.character(network_mod[[1]]),
                                         as.character(network_mod[[2]])),
                              modularity = as.character(network_mod[[3]]))
    # Add category of the node
    network_lab$cat <- ifelse(network_lab$node %in% rownames(dat),
                              "site", "sp")

  } else if(algo %in% c("greedy", "girvan", "walktrap", "louvain", "infomap")){
    # https://stats.stackexchange.com/questions/209086/community-detection-and-modularity
    # https://www.sixhat.net/finding-communities-in-networks-with-r-and-igraph.html

    # Transforming matrix into square matrix (first pixels and then species)
    rownames_mat <- rownames(dat)
    dat_sq <- rbind(
      cbind(array(0, c(nrow(dat), nrow(dat))), as.matrix(dat)),
      cbind(as.matrix(t(dat)), array(0, c(ncol(dat), ncol(dat)))))

    # Add pixel names to square matrix
    rownames(dat_sq)[1:length(rownames_mat)] <- rownames_mat
    colnames(dat_sq)[1:length(rownames_mat)] <- rownames_mat

    if(weight == FALSE){
      # Tranforming site_sp matrix into binary matrix
      dat_sq[dat_sq > 0] <- 1
      network <- graph_from_adjacency_matrix(
        dat_sq, mode = "undirected", add.rownames = NULL, weighted = NULL)
    } else if(weight == TRUE){
      network <- graph_from_adjacency_matrix(
        dat_sq, mode = "undirected", add.rownames = NULL, weighted = TRUE)
    }

    # Modularity algorithm
    if(algo == "greedy"){
      network_mod <- fastgreedy.community(network)
    } else if(algo == "girvan"){
      network_mod <- cluster_edge_betweenness(network, modularity = TRUE)
    } else if(algo == "walktrap"){
      network_mod <- cluster_walktrap(graph = network)
    } else if(algo == "louvain"){
      network_mod <- cluster_louvain(graph = network)
    } else if(algo == "infomap"){
      network_mod <- cluster_infomap(graph = network)
    }
    # Convert results into data.frame
    network_lab <- c()
    for(i in 1:length(igraph::groups(network_mod))){
      network_lab <- rbind(network_lab,
                           cbind(igraph::groups(network_mod)[[i]],
                                 rep(as.numeric(i),
                                     length(igraph::groups(network_mod)[[i]]))))
    }
    network_lab <- data.frame(network_lab)
    colnames(network_lab) <- c("node", "module")
    # rownames(network_lab) <- network_lab$node

    # Add modularity score
    network_lab$modularity <- modularity(network_mod)

    # Add category of the node
    network_lab$cat <- ifelse(network_lab$node %in% rownames(dat),
                              "site", "sp")
  }

  return(network_lab)

  # Also try with bipartite package
  # computeModules(web, method = "Beckett", deep = FALSE,
  #                deleteOriginalFiles = TRUE, steps = 1000000, tolerance = 1e-10,
  #                experimental = FALSE, forceLPA = FALSE)

  # And rnetcarto
}
