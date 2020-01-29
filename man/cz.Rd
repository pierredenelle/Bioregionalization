\name{cz}
\alias{cz}
\title{Calculate the Cz scores of sites and species for a given
bipartite network.}
\usage{
cz(dat, sp_col, site_col, bip, ab = NULL)
}
\arguments{
\item{dat}{a data.frame with each row indicating the presence
of a species in a given site.}

\item{bip}{bip contains three columns: node, mod and cat which respectively
stand for the name of the node, its module and its bipartite category.}

\item{sp_col}{a character string corresponding to the column
    with species.}

\item{site_col}{a character string corresponding to the column
    with sites.}

\item{ab}{a logical indicating whether abundances should be taken into
account.}

}
\description{
The function returns a list containing two data.frame which correspond to
bip data.frame with extra-columns Cz, z being only standardized per species
in the second data.frame of the list.

See Guimera et Amaral, Nature, 2005.
}
\details{
# column of abundance can be added to weight the calculation
# (not implemented yet)
}
\section{Note}{}

\examples{
}
