\name{comparison}
\alias{comparison}
\title{Function to identify pixels grouped together through different
bioregionalizations.}
\usage{
all_maps(dat, site, sp, ab, binary = FALSE)
}
\arguments{
\item{dat}{a data.frame with a column containing the sites and
columns containing different partitions.}

\item{bio_col}{a numeric indicating which columns of dat contains
the different bioregions.}

\item{thres}{a numeric determining the minimum number of pixels per pair.}

}
\description{
Function to identify pixels grouped together through different methods of
bioregionalizations.

}
\details{}
\section{Note}{}

\examples{
ex <- data.frame(site = paste0("site", seq(1:7)),
bio1 = c("a", "b", "a", "c", "d", "b", "d"),
bio2 = c("aa", "aa", "aa", "bb", "cc", "bb", "cc"))

comparison(ex, bio_col = c(2, 3), thres = 1)
}
