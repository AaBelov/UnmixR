##' unmixR implements N-FINDR and Vertex Component Analysis (VCA)
##' algorithms which can recover pure component spectra and their
##' respective concentrations from a hyperspectral data set.
##'
##' @name unmixR-package
##' @title Hyperspectral Unmixing Methods
##' @docType package
##'
##' @author Conor McManus
##'
##' Maintainer: Claudia Beleites <chemometrie@beleites.de>
##' @rdname unmixR-package
##' @aliases unmixR
##' @keywords package
##' @import hyperSpec
##' @import nnls
##' @import MASS
# @import svUnit

# This code needs to be here due to the order in which
# things are sourced, apparently.  It made more sense to
# me to put it unittests.R but that causes problems.  BH.
	
if (!requireNamespace ("svUnit", quietly = TRUE)){
  `.test<-` <- function (f, value) {
    class (value) <-  c ("svTest", "function")
    attr (f, "test") <- value
    f
  }
} else {
  `.test<-` <- svUnit::`test<-`
}
