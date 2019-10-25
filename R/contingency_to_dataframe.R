
contingency_to_dataframe <- function(dat){
  require(DescTools)

  if(!is.matrix(dat)){
    stop("Contingency table should be a matrix with sites as rows and
         species as columns.")
  }

  dat_df <- as.data.frame(as.table(dat))
  colnames(dat_df) <- c("site", "sp", "ab")
  dat_df <- dat_df[which(dat_df$ab > 0), ]

  return(mat)
}
