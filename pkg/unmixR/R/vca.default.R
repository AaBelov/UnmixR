##' @name vca
##' @rdname vca
##' @export
##' @include unmixR-package.R
##' @include vca.R

vca.default <- function(data, p, method= c("Mod", "Lopez", "05"), seed=NULL, ...) {

  # check if the method passed in is valid
  method <- match.arg (method)

  # transform the input into a matrix
  data <- as.matrix (data)

  # check for p being with the valid range, >= 2
  if (!is.numeric (p) || p < 2 || p > ncol (data)) {
    stop("p must be a positive integer >= 2 and <= ncol (data)")
  }

  # set the random number generator seed if supplied
  if (!is.null(seed)) {
    set.seed(seed)
  }

  vcaFunc <- get(paste("vca", method, sep=""))
  val <- vcaFunc(data, p, ...)

  struct <- list(
    data = data,
    indices = sort(val) # TODO: sort should not be necessary at this point 
    # as all functions already sort 
    )

    structure(struct, class = "vca")
}


.test(vca.default) <- function() {
  # Note: .testdata$x matches all columns of .testdata, which are x.L1, x.L2, x.L3
  # test: vca produces error for invalid values of p
  checkException (vca (.testdata$x, p = "---"))
  checkException (vca (.testdata$x, p = 0))
  checkException (vca (.testdata$x, p = 1))
  checkException (vca (.testdata$x, p = 4))

  # test: vca produces error for invalid method
  checkException(vca (.testdata$x, p, method="invalid"))

  # test correct calculations for the available methods

  methods <- eval (formals (vca.default)$method)

  for (m in methods) {
    ## .testdata has 3 components, so picking 2 out of 3
    model <- vca (.testdata$x, p = 2, method = m)
    checkTrue (all (model$indices %in% .correct),
               msg = sprintf ("%s: .testdata, p = 2 yields %s", m, paste (model$indices, collapse = " ")))
    
    # BH changed to p = 3, since once in VCA05 this data passes to the 2nd block where d = p - 1 
    checkTrue (all (vca (.testdata$x, p = 3, method = m)$indices %in% .correct), 
               msg = sprintf ("%s: .testdata, p = 2", m))
    
    ## all 3 components should be recovered, vca output is sorted.
    model <- vca (.testdata$x, p = 3, method = m)
    checkEquals (model$indices, .correct,
                 msg = sprintf ("%s: .testdata, p = 2 yields %s", m, paste (model$indices, collapse = " ")))
  }

  # test: if hyperSpec is available, test on hyperSpec object
  # tests also the correct application of as.matrix.
  if (require ("hyperSpec")) {
    checkEquals (vca (laser, p = 2)$indices, c (4, 81))
  }
}
