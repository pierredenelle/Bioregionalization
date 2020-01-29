\name{zscore}
\alias{zscore}
\title{Compute zscores for species from a given delineation of sites.}
\usage{
zscore(dat, sp_col, site_col, bioregion_col)
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

\item{output_format}{a character string equal to 'matrix' or 'data.frame'
that determines the format of the output}

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
                   bioregion_col = "bioregion")
z_scores
}
