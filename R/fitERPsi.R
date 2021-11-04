#' @include entryFit.R unblindFit.R psiFit.R
.fitERPsi <- function(data, 
                      entry,  
                      unblind0,  
                      unblind1,  
                      Psi,  
                      dObj,
                      infTimesBlinded,
                      wgtType) {

  if (wgtType > 0L) message("performing regressions")

  if (wgtType == 2L) {

    # entryFit returns the vector {n} f_E
    fEX <- .entryFit(data = data, model = entry, dObj = dObj)

    # unblindFit returns a list containing:
    #   "unblindFit0" the returned fit object for the A = 0 analysis
    #   "unblindFit1" the returned fit object for the A = 1 analysis
    #   "fR" a matrix {n x 2} of f_R for A = 0/1
    #   "LambdaT" a matrix {nTimesBlinded x 2} Lambda(u) for A = 0/1
    #   "expXbeta" a matrix {n+2 x 2} exp(x beta) with that for the mean for
    #              A = 0 and the mean for A = 1 tacked on to the end.
    unblindFit <- .unblindFit(data = data, 
                              modelA0 = unblind0,  
                              modelA1 = unblind1,  
                              dObj = dObj,  
                              times = infTimesBlinded)

    # psiFit returns the vector {n} p_{Psi}
    Psiprob <- .psiFit(data = data, model = Psi, dObj = dObj)

  } else {
    # Rcpp does not like NULLs
    fEX <- 0.0
    unblindFit <- list("fR" = matrix(data = 0.0, nrow = 1L, ncol = 1L), 
                       "LambdaT" = matrix(data = 0.0, nrow = 1L, ncol = 1L), 
                       "expXbeta" = matrix(data = 0.0, nrow = 1L, ncol = 1L))
    Psiprob <- 0.0
  }

  return( list("fEX" = fEX, "unblindFit" = unblindFit, "Psiprob" = Psiprob) )
}
