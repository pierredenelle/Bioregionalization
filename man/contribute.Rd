\name{contribute}
\alias{contribute}
\title{Compute Cz and rho scores for species from a given delineation of
sites.}
\usage{
contribute(dat, sp_col, site_col, bioregion_col)
}
\arguments{
\item{dat}{a data.frame with each row indicating the presence
of a species in a given site.}

\item{sp_col}{a character string corresponding to the column
    with species.}

\item{site_col}{a character string corresponding to the column
    with sites.}

\item{bioregion_col}{a character string corresponding to the column
    with bioregions.}

\item{bioregion_sp_col}{a character string corresponding to the column
    with bioregions of species. Can be null in case of a projected
    network.}

\item{ab}{a logical indicating whether abundances should be taken into
account.}

}

\description{
See Guimera et Amaral, Nature, 2005.
See Lenormand et al., Ecology and Evolution, 2019.
}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
species = c("a", "b", "a", "c", "d", "b", "d"),
count = c(10, 100, 1, 20, 50, 10, 20),
bioregion = c(rep("region1", 5), rep("region2", 2)))

scores <- contribute(dat = ex, sp_col = "species", site_col = "sites",
                 bioregion_col = "bioregion")
scores

}
