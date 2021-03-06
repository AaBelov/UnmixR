##' Michael E. Winter's 1999 N-FINDR unmixing algorithm
##' 
##' This technique is based on the fact that, in N spectral dimensions, the
##' N-volume contained by a simplex formed of the purest pixels is larger
##' than any other volume formed from any other combination of pixels.
##' Intended to be called from \code{\link{nfindr}}.
##' 
##' @param data Data matrix to unmix
##' @param p Number of endmembers
##' @param indices Indices used in the simplex estimation
##' @param iters Max number of iterations, defaults to 3*p
##' @include unmixR-package.R
##' @return The indices that indicate the position of the endmembers in the
##'   original dataset
##' @export
##'
##' @references Michael E. Winter; "N-FINDR: an algorithm for fast autonomous
##'   spectral end-member determination in hyperspectral data", Proc.
##'   SPIE 3753, Imaging Spectrometry V, 266 (October 27, 1999), 
##'   doi:10.1117/12.366289
##' @export
nfindr99 <- function(data, p, indices, iters=3*p) {
  simplex <- .simplex(data, p, indices)
  nspectra <- nrow(data)

  # calculate the initial volume using the random endmembers
  volume <- abs(det(simplex))
  volume.prev <- -1
  volume.now <- volume

  # keep replacing endmembers until there is never an increase in volume
  # or the max iterations are reached (indicates pure endmembers not found)
  iter <- 1
  while (volume.now > volume.prev && iter <= iters) {
    for (k in 1:p) {
      for (i in 1:nspectra) {
        # store current sample as it may need to be placed back into the
        # simplex after the following replacement
        sample <- simplex[2:p,k]

        # replace the k-th endmember with the i-th reduced spectrum
        # and recalculate the volume
        simplex[2:p,k] <- data[i,]
        testVolume <- abs(det(simplex))

        # if the replacement increased the volume then keep the replacement
        # and the note the spectrum's index
        if (testVolume > volume) {
          volume <- testVolume
          indices[k] <- i
        }
        # otherwise revert the replacement
        else {
          simplex[2:p,k] <- sample
        }
      }
    } # end of (k in 1:p) loop
    
    iter <- iter+1
    volume.prev <- volume.now
    volume.now <- volume
  } # end of 'while'
  
  indices <- sort(indices)
  indices
}
