
contingency_to_dataframe <- function(dat, col1 = "site", col2 = "sp",
                                     col3 = "ab", remove_zeros = TRUE){
  require(DescTools)
  
  if(!is.matrix(dat)){
    stop("Contingency table should be a matrix with sites as rows and
         species as columns.")
  }
  
  if(!is.character(col1)){
    stop("col1 must be the column name describing the rows of the contingency
         matrix.")
  }
  
  if(!is.character(col2)){
    stop("col2 must be the column name describing the columns of the
    contingency matrix.")
  }
  
  if(!is.character(col3)){
    stop("col3 must be the column name describing the values of the contingency
         matrix.")
  }
  
  if(!is.logical(remove_zeros)){
    stop("remove_zeros must be a boolean determining whether the null values
         from the contingency matrix has to be removed from the output.")
  }
  
  # Conversion as data.frame
  dat_df <- as.data.frame(as.table(dat))
  colnames(dat_df) <- c(col1, col2, col3)
  
  if(remove_zeros == TRUE){
    dat_df <- dat_df[which(dat_df$ab > 0), ]
  }
  
  return(dat_df)
}
