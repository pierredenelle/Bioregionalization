\name{project_network}
\alias{project_network}
\title{Project a contingency table using a similarity metric}
\usage{
project_network(contingency_mat, similarity = "simpson")
}
\arguments{
\item{contingency_mat}{a contingency table with sites as rows
and species as columns}

\item{similarity}{a character string indicating which similarity metric
to chose among simpson, jaccard, sorensen and whittaker}

}
\description{
Project a contingency table into a similarity data.frame,
the similarity being based upon the sites each pair of species has in common.
}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3)),
species = c("a", "b", "a", "c", "d"),
count = c(10, 100, 1, 20, 50))

ex_mat <- contigency(dat = ex, sp = "species", site = "sites", ab = "count")
project_network(ex_mat, similarity = "simpson")
}
