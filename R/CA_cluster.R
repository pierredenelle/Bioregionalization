
CA_cluster <- function(dat){
  require(ade4)

  # CA analysis
  CA_site_sp <- dudi.coa(dat, scannf = FALSE, nf = 2)
  # plot(CA_site_sp$li, pch = 16, col = rgb(0, 0, 0, 0.5))

  # K-means clustering
  # Determine number of clusters
  wss <- (nrow(dat) - 1) * sum(apply(dat, 2, var))
  for (i in 2:15){
    wss[i] <- sum(kmeans(dat, centers = i)$withinss)
  }
  # plot(1:15, wss, type = "b", xlab = "Number of Clusters",
  #      ylab = "Within groups sum of squares")

  # Number of cluster centers determined after evaluation of elbow plots
  k_mean <- kmeans(CA_site_sp$li, centers = 17) # number of centers?
  #plot(CA_site_sp$li, col = po$clust, pch = 16)

  res <- data.frame(site = names(k_mean$cluster),
                    cluster = as.character(k_mean$cluster))

  return(res)
}

# With FactoMineR
# library(FactoMineR)
# CA_site_sp2 <- CA(dat, graph = FALSE)
# CA_hcpc <- HCPC(CA_site_sp2, nb.clust = 3, graph = FALSE)
# plot(CA_hcpc, axes = 1:2)
#
# str(CA_hcpc_hcpc$desc.ind)
# head(CA_hcpc_hcpc$desc.var$`1`) # cluster for species
# CA_hcpc_hcpc$desc.ind # cluster for communities
