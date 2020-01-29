\name{interact}
\alias{interact}
\title{Computes interactions between modules.}
\usage{
interact(dat, plot = FALSE, output_format = "matrix")
}
\arguments{
\item{input_network}{a character string equal to 'projected' or
'bipartite'. It determines the format of the input.}

\item{dat}{when input_network is 'projected', a data frame containg species in
rows, bioregions in columns and filled with the zscore of each species/
bioregion combination. Can be generated with the function zscore.
When input_network is 'bipartite', dat must be a data.frame with four columns
containing the links between the two categories of nodes and their respective
modules.}

\item{plot}{a character string corresponding to the column
    with sites.}

\item{output_format}{a boolean determining whether the barplot of interactions
between bioregions should be saved.}

}
\description{
See Lenormand et al., Ecology and Evolution, 2019.
}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20),
bioregion = c(rep("region1", 5), rep("region2", 2)))

z_scores <- zscore(ex, sp_col = "species", site_col = "sites",
                   bioregion_col = "bioregion", output_format = "matrix")
lambda <- interact(input_network = "projected",
dat = z_scores, plot = TRUE, output_format = "matrix")
lambda
}
