#'
#' Virtual dataset used for the vignette
#'
#' Random raster with degree of autocorrelation generated
#' 100 virtual species with different response curve to the gradient
#' Mean and standard deviation of species vary, such as some of them are more
#' or less generalist/specialist.
#' From these species, suitability index is drawn and is used as a proxy of
#' local abundance.
#'
#' @docType data
#'
#' @usage data(virtual_sp)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#' @references Leroy et al. (2016) Ecography 39: 599-607
#' (\href{https://onlinelibrary.wiley.com/action/showCitFormats?doi=10.1111%2Fecog.01388}{PubMed})
#'
#' @source \href{https://phenome.jax.org/projects/Moore1b}{QTL Archive}
#'
#' @examples
#' data(virtual_sp)
"virtual_sp"
