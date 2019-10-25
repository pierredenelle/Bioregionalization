\name{contingency_to_dataframe}
\alias{contingency_to_dataframe}
\title{Create a data.frame from a contingency table.}
\usage{
contingency_to_dataframe(dat)
}
\arguments{
\item{dat}{a contingency matrix with sites as rows and species as columns.}

}
\description{
Converts a contingency table into a data.frame with replicated rows per site.
}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20))

mat <- contingency(dat = ex, site = "sites", sp = "species", ab = "count")
contingency_to_dataframe(mat)
}
