
ward_cluster <- function(dat, method = "ward.D2", optim_method = "firstSEmax",
                         nstart = 25, B = 50, K.max = 20){
  require(cluster)

  if(!is.matrix(dat)){
    stop("dat should be a matrix with sites as rows and species as columns.")
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

  if(K.max > nrow(dat)){
    stop("K.max should not be superior to the number of rows of the
         contingency matrix.")
  }

  if(K.max > nrow(unique(dat))){
    stop("K.max should not be superior to the unique number of rows of the
         contingency matrix.")
  }

  # Euclidean distance matrix
  # dist_sp_mat <- dist(dat)
  euc_dist <- function(m) {
    mtm <- Matrix::tcrossprod(m)
    sq <- rowSums(m*m)
    sqrt(outer(sq,sq,"+") - 2*mtm)
  } 
  
  dist_sp_mat <- euc_dist(dat)
  dist_sp_mat <- as.dist(dist_sp_mat)
  
  #h <- hclust(dist_sp_mat, method = method)
  require(fastcluster)
  h <- fastcluster::hclust(dist_sp_mat, method = method)
  # plot(h)

  # Determine optimal numbers of clusters
  # https://stackoverflow.com/questions/53159033/how-to-get-the-optimal-number-of-clusters-from-the-clusgap-function-as-an-output
  # https://uc-r.github.io/kmeans_clustering#gap
  # https://uc-r.github.io/hc_clustering

  gap_stat <- clusGap(dat, FUN = kmeans, nstart = nstart, K.max = K.max,
                      B = B)

  optim_k <- maxSE(f = gap_stat$Tab[, "gap"],
                   SE.f = gap_stat$Tab[, "SE.sim"],
                   method = optim_method)

  # Cut the tree with optim_k numbers
  dend <- as.dendrogram(h)
  # Data.frame of results
  res <- data.frame(site = names(dendextend::cutree(dend,
                                                    k = optim_k)),
                    cluster = as.character(dendextend::cutree(dend,
                                                              k = optim_k)))

  return(res)
  # Visualization
  # factoextra::fviz_nbclust(dat, kmeans, method = "gap_stat", k.max = 20)

  # dend <- as.dendrogram(h)
  # library(dendextend)
  # labels(dend) <- x[order.dendrogram(dend)]
  # # due to the ties - there is specific reason to have this be these 3 clusters:
  # cutree(dend, k = 3)
  #
  # library(NbClust)
  # NbClust(data = dat, diss = NULL, distance = "euclidean",
  #         min.nc = 1, max.nc = 15, method = "ward.D2")
}

