\name{all_maps}
\alias{all_maps}
\title{Wrap-up function that runs several methods available in the package.}
\usage{
all_maps(dat, site, sp, ab, binary = FALSE)
}
\arguments{
\item{dat}{a data.frame with duplicated rows for sites, one row per species,
a third column indicating the abundance of the species in this site}

\item{site}{the character string corresponding to the site column name}

\item{sp}{the character string corresponding to the species column name}

\item{ab}{the character string corresponding to the abundance column name}

\item{binary}{a logical indicating if abundances should be converted as 0/1}

}
\description{
Wrap-up function that runs several methods available in the package.
}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20))


}
