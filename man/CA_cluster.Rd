\name{CA_cluster}
\alias{CA_cluster}
\title{Performs a CA analysis and k-means clustering from a contingency table}
\usage{
contingency(dat, sp, site, ab, binary = FALSE)
}
\arguments{
\item{dat}{a contingency table}

}
\description{
Performs a CA analysis and k-means clustering from a contingency table and
extract clusters.
}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20))

tmp <- contingency(dat = ex, sp = "species", site = "sites", ab = "count")
CA_cluster(tmp)
}
