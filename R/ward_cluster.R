
ward_cluster <- function(dat){
  require(cluster)
  # Distance matrix
  dist_sp_mat <- dist(dat)
  h <- hclust(dist_sp_mat, method = "ward.D2")
  # plot(h)

  # Determine optimal numbers of clusters
  # https://stackoverflow.com/questions/53159033/how-to-get-the-optimal-number-of-clusters-from-the-clusgap-function-as-an-output
  # https://uc-r.github.io/kmeans_clustering#gap
  # https://uc-r.github.io/hc_clustering

  gap_stat <- clusGap(sp_mat, FUN = kmeans, nstart = 25,
                      K.max = 20, B = 50)

  optim_k <- maxSE(f = gap_stat$Tab[, "gap"],
                   SE.f = gap_stat$Tab[, "SE.sim"],
                   method = "firstSEmax")

  # Cut the tree with optim_k numbers
  dend <- as.dendrogram(h)
 # Data.frame of results
  res <- data.frame(site = names(dendextend::cutree(dend,
                                                    k = optim_k)),
                    cluster = as.character(dendextend::cutree(dend,
                                                              k = optim_k)))

  return(res)
  # Visualization
  # factoextra::fviz_nbclust(sp_mat, kmeans, method = "gap_stat", k.max = 20)

  # dend <- as.dendrogram(h)
  # library(dendextend)
  # labels(dend) <- x[order.dendrogram(dend)]
  # # due to the ties - there is specific reason to have this be these 3 clusters:
  # cutree(dend, k = 3)
  #
  # library(NbClust)
  # NbClust(data = sp_mat, diss = NULL, distance = "euclidean",
  #         min.nc = 1, max.nc = 15, method = "ward.D2")
}

