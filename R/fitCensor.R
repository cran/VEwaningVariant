#' @include censorFit.R
.fitCensor <- function(data, 
                       censor0, 
                       censor1, 
                       dObj, 
                       wgtType, 
                       infTimesBlinded, 
                       infTimesUnblinded, 
                       L) {

  if (wgtType > 0L) {

    # CensorFit returns a list containing
    #    "dropFit0" the regression object returned for the A = 0 analysis
    #    "dropFit1" the regression object returned for the A = 1 analysis
    #    "LambdaTU" {nTimesU x 2} Lamdba(u) on timesU for A = 0/1
    #    "LambdaTB" {nTimesB x 2} Lamdba(u) on timesB for A = 0/1
    #    "LambdaR" {n+2 x 2} Lambda(r) for A = 0/1 with means tacked on
    #    "expXbetaG0" {n+2 x 2} exp(xbeta) when G = 0 for A = 0/1
    #    "expXbetaG1" {n+2 x 2} exp(xbeta) when G = 1 for A = 0/1
    censorFit <- .censorFit(modelA0 = censor0, 
                            modelA1 = censor1,  
                            data = data,  
                            dObj = dObj,  
                            timesB = infTimesBlinded, 
                            timesU = infTimesUnblinded,
                            L = L)

  } else {
    censorFit <- list("LambdaTU" = matrix(data = 0.0, nrow = 1L, ncol = 1L), 
                      "LambdaTB" = matrix(data = 0.0, nrow = 1L, ncol = 1L), 
                      "LambdaR" = matrix(data = 0.0, nrow = 1L, ncol = 1L),
                      "expXbetaG0" = matrix(data = 0.0, nrow = 1L, ncol = 1L),
                      "expXbetaG1" = matrix(data = 0.0, nrow = 1L, ncol = 1L))
  }

  return( censorFit )
}
