######################################################################

##  Fit Cox models for unblinding hazards for Gam = 1 and A = 0, A = 1
##  calculate survival probabilities K_R at U's for all individuals
##  with Delta = 1 whose U's fall in [TP,TU), Gam=1 and [TU,TC), Gam=2.
##  Also calculate the survival probabilities at the mean X values to
##  form stabilized weights and get the f_R at the R values

##  merge() used below will order the KR1 and KR2 data frames by ID
##  number in original data set
#
#' @importFrom stats model.response
#' @importFrom stats model.frame
#' @import survival
#'
# @param data A data.frame object. All time-independent data required for the 
#   analysis.
#
# @param tdData A data.frame object or NULL. All time-dependent data.
#
# @param modelA0 A formula object. The coxph model for R with Gamma = 1
#   and A = 0. Note that the LHS is taken as a Surv() object as defined by 
#   package survival.
#
# @param modelA1 A formula object. The coxph model for R with Gamma = 1
#   and A = 1. Note that the LHS is taken as a Surv() object as defined by 
#   package survival.
#
# @param dObj A list object. The column headers of data containing required
#   variables
#
# @param times A numeric vector. The times at which Lambda is to be calculated
#
# returns a list containing:
#   "unblindFit0" the returned fit object for the A = 0 analysis
#   "unblindFit1" the returned fit object for the A = 1 analysis
#   "fR" a matrix {n x 2} of f_r for A = 0/1
#   "LambdaT" a matrix {nTimesBlinded x 2} Lambda(u) for A = 0/1
#   "expXbeta" a matrix {n+2 x 2} exp(x beta) with that for the mean for
#              A = 0 and the mean for A = 1 tacked on to the end.
#
#' @import survival
#' @importFrom stats reformulate predict coef terms
#
#' @include survFunc.R pzero.R
.unblindFit <- function(data, 
                        modelA0, 
                        modelA1,
                        dObj,
                        times) {

  message("\tunblinding")

  n <- nrow(x = data)

  # update formulae to include Surv() LHS
  survForm <- paste0('Surv(',dObj$R, ",", dObj$Gam,')')

  modelA0 <- stats::reformulate(termlabels = attr(x = terms(x = modelA0), 
                                                  which = "term.labels"),
                                response = survForm, 
                                intercept = FALSE)

  modelA1 <- stats::reformulate(termlabels = attr(x = terms(x = modelA1), 
                                                  which = "term.labels"),
                                response = survForm, 
                                intercept = FALSE)

  #  Fit Cox model for R using only A = 0 participants
  subsetA0 <- data[,dObj$A] == 0L

  A0Fit <- tryCatch(expr = survival::coxph(formula = modelA0, 
                                           data = data,
                                           subset = subsetA0),
                    error = function(e) {
                              stop("unable to fit Cox model for unblinding, A = 0\n",
                                   e$message, call. = FALSE)
                            })

  if (anyNA(x = coef(object = A0Fit))) {
    stop("fit of Cox model for unblinding, A = 0 resulted in NA coefficients",
         call. = FALSE)
  }

  #  Fit Cox model for R using only A = 1 participants
  subsetA1 <- data[,dObj$A] == 1L
  A1Fit <- tryCatch(expr = survival::coxph(formula = modelA1, 
                                           data = data,
                                           subset = subsetA1),
                    error = function(e) {
                              stop("unable to fit Cox model for unblinding, A = 1\n",
                                   e$message, call. = FALSE)
                            })

  if (anyNA(x = coef(object = A1Fit))) {
    stop("fit of Cox model for unblinding, A = 1 resulted in NA coefficients",
         call. = FALSE)
  }

  ## Predictions

  ##  Form stabilized f_R densities evaluated at the R values

  KRr.0 <- .survFunc(object = A0Fit, newdata = data) * 
           .pzero(object = A0Fit, newdata = data)

  KRr.1 <- .survFunc(object = A1Fit, newdata = data) * 
           .pzero(object = A1Fit, newdata = data)

  ##  For stabilized weights must get the same thing for the means of X

  meansA0 <- colMeans(x = data[subsetA0,,drop=FALSE])

  newdata <- data
  for (i in 1L:length(x = meansA0)) {
    newdata[,i] <- meansA0[i]
  }
  newdata[,dObj$R] <- data[,dObj$R]
  newdata[,dObj$Gam] <- data[,dObj$Gam]

  KRr.mean.0 <- .survFunc(object = A0Fit, newdata = newdata) * 
                .pzero(object = A0Fit, newdata = as.data.frame(x = t(x = meansA0)))

  meansA1 <- colMeans(x = data[subsetA1,,drop=FALSE])

  for (i in 1L:length(x = meansA1)) {
    newdata[,i] <- meansA1[i]
  }
  newdata[,dObj$R] <- data[,dObj$R]
  newdata[,dObj$Gam] <- data[,dObj$Gam]

  KRr.mean.1 <- .survFunc(object = A1Fit, newdata = newdata) * 
                .pzero(object = A1Fit, newdata = as.data.frame(x = t(x = meansA1)))

  ## need Lamba(t) and exp(Xbeta) for each treatment
  # ts are the infection times during the blinded stage
  # exp(Xbeta) includes the average values for each treatment group

  newdata <- data
  newdata <- rbind(newdata, meansA0, meansA1)

  if (!is.null(x = times)) {
    LambdaT0 <- suppressWarnings(expr = summary(object = survfit(formula = A0Fit),
                                                time = times,
                                                extend = TRUE)$cumhaz)
    LambdaT1 <- suppressWarnings(expr = summary(object = survfit(formula = A1Fit),
                                                time = times,
                                                extend = TRUE)$cumhaz)
  } else {
    LambdaT0 <- 0.0
    LambdaT1 <- 0.0
  }

  expXbeta0 <- predict(object = A0Fit, newdata = newdata, type = "risk")

  expXbeta1 <- predict(object = A1Fit, newdata = newdata, type = "risk")

  return( list("unblindFit0" = A0Fit,
               "unblindFit1" = A1Fit,
               "fR" = cbind(KRr.mean.0 / KRr.0, KRr.mean.1 / KRr.1),
               "LambdaT" = cbind(LambdaT0, LambdaT1),
               "expXbeta" = cbind(expXbeta0, expXbeta1)) )
}
