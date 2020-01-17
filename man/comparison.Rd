\name{comparison}
\alias{comparison}
\title{Function to identify pixels grouped together through different bioregionalizations.}
\usage{
all_maps(dat, site, sp, ab, binary = FALSE)
}
\arguments{
\item{dat}{a data.frame with duplicated rows for sites, one row per species,
a third column indicating the abundance of the species in this site}

\item{site}{the character string corresponding to the site column name}

}
\description{
Function to identify pixels grouped together through different bioregionalizations.
}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20))


}
