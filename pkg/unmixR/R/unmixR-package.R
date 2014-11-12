##' Implements the N-FINDR and Vertex Component Analysis algorithms
##' to recover pure component spectra and their respective concentrations
##' from a set of measured spectra.
##' 
##' @name unmixR-package
##' @title Spectral Unmixing Methods
##' @docType package
##' @author Conor McManus
##' 
##' Maintainer: Conor McManus <chathurga@@gmail.com>
##' @rdname unmixR-package
##' @keywords package
{
  if (!require ("svUnit", quietly = TRUE)){
    `.test<-` <- function (f, value) {
      class (value) <-  c ("svTest", "function")
      attr (f, "test") <- value
      f
    }
  } else {
    `.test<-` <- svUnit::`test<-`
  }
}