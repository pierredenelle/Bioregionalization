
project_network <- function(contingency_mat, similarity = "simpson"){
  require(Rcpp)
  require(SMUT)

  if(!is.matrix(contingency_mat)){
    stop("Contingency table should be a matrix with sites as rows and
         species as columns.")
  }

  if(!(similarity %in% c("simpson", "jaccard", "sorensen", "whittaker",
                         "bray"))){
    stop("Similarity metric chosen is not available.
     Please chose among the followings:
         simpson, jaccard, sorensen, whittaker or bray")
  }

  # Convert as matrix to remove names
  contingency_mat <- as.matrix(contingency_mat)

  if(similarity == "bray"){
    require(ecodist)
    require(reshape2)
    bray <- bcdist(contingency_mat, rmzero = FALSE)
    # Conversion as matrix and removal of upper part and diagonal
    bray <- as.matrix(bray)
    bray[upper.tri(bray)] <- NA
    diag(bray) <- NA

    # Convert distance matrix into data.frame
    abc <- melt(bray, varnames = c("id1", "id2"))
    # Remove NAs
    abc <- abc[complete.cases(abc), ]
    colnames(abc) <- c("id1", "id2", "bray")
    # Create columns of site names and put id in the id columns
    if(!is.null(rownames(contingency_mat))){
      abc$id1_name <- abc$id1
      abc$id2_name <- abc$id2
      abc$id1 <- match(abc$id1, rownames(contingency_mat))
      abc$id2 <- match(abc$id2, rownames(contingency_mat))
    }
    # Similarity instead of dissimilarity
    abc$bray <- 1 - abc$bray

    # Remove zeros
    abc <- abc[which(abc$bray > 0), ]

  } else{

    # Get the number of species shared by two sites
    # eigenMapMatMult to compute matrix product
    # equivalent to contingency_mat %*% t(contingency_mat)
    a <- eigenMapMatMult(contingency_mat, t(contingency_mat))

    abc <- which(a > 0, arr.ind = TRUE)
    abc <- cbind(abc, a[a > 0])
    # Remove upper part of the matrix
    abc <- abc[abc[, 1] < abc[, 2], ]

    a <- diag(a) # number of species per sites
    b <- a[abc[, 1]] # number of species in first column
    c <- a[abc[, 2]] # number of species in second column

    # Bind informations
    abc <- cbind(abc, b, c)
    # b equals to species only in first column
    abc[, 4] <- abc[, 4] - abc[, 3]
    # c equals to species only in second column
    abc[, 5] <- abc[, 5] - abc[, 3]

    colnames(abc) <- c("id1", "id2", "a", "b", "c")

    # Conversion as data.frame
    abc <- as.data.frame(abc)

    # Similarity metric
    if(similarity == "simpson"){
      abc$simpson <- 1 - pmin(abc$b, abc$c)/(abc$a + pmin(abc$b, abc$c))
    } else if(similarity == "sorensen"){
      abc$sorensen <- 2*abc$a/(2*abc$a + abc$b + abc$c)
    } else if(similarity == "jaccard"){
      abc$jaccard <- abc$a/(abc$a + abc$b + abc$c)
    } else if(similarity == "whittaker"){
      abc$whittaker <- (abc$a + abc$b + abc$c)/((2*abc$a + abc$b + abc$c)/2)
    }
  }
  # If contingency matrix has rownames, ressign them to abc data.frame
  if(!is.null(rownames(contingency_mat))){
    abc$id1_name <- rownames(contingency_mat)[abc$id1]
    abc$id2_name <- rownames(contingency_mat)[abc$id2]
  }

  return(abc)
}
