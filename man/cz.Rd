\name{run_OSLOM}
\alias{run_OSLOM}
\title{Convert the .tp output of OSLOM algorithm into a data.frame,
with each site having its bioregion.}
\usage{
run_OSLOM(dat, n_runs = 1, t_param = 0.1, cp_param = 0.5, saving_directory)
}
\arguments{
\item{dat}{a .tp object generated by OSLOM algorithm. Import in R with
readLines() function. This object can be generated with run_oslom() function.}

\item{n_runs}{Numeric parameter setting the number of runs for the first
hierarchical level. The default value is 10. The bigger, the slower the
function but the more accurate the results.}

\item{t_param}{Numeric parameter setting the probability value to consider
a module as significant. Must be comprised between 0 and 1.}

\item{cp_param}{Numeric parameter setting the resolution. Must be comprised
between 0 and 1.}

\item{saving_directory}{Path where the OSLOM .tp file containing
         the bioregions identified will be saved.}

}
\description{
Run OSLOM algorithm and save .tp file with the identified bioregions in a
chosen directory.
See http://www.oslom.org/code/ReadMe.pdf for more precisions.
}
\details{}
\section{Note}{}

\examples{
}