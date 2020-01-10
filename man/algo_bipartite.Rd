\name{algo_bipartite}
\alias{algo_bipartite}
\title{Bipartite community detection}
\usage{
algo_bipartite(dat, algo = "greedy", weight = FALSE)
}
\arguments{
\item{dat}{a contingency matrix with sites as rows and species as columns.}

\item{algo}{An algorithm used to perfom the bipartite community detection.
Available options are greedy, girvan, walktrap, louvain, infomap, beckett and dorman.}

\item{weight}{A boolean indicating whether the detection of communities has to
be weighted or not.}

}
\description{
Run a community detection algorithm on a bipartite network.
Algorithms available are from igraph package (greedy, girvan, walktrap, louvain, infomap)
and also from the bipartite package (beckett and dorman).
Algorithm can be weighted in case of a weighted bipartite network.
}
\details{}
\section{Note}{}

\examples{
}
