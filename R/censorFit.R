######################################################################
##  Fit Cox models for Censor hazards
##  merge() used below will order the KR data frames by ID
##  number in original data set
#
# @param modelA0 A formula object. The coxph model for Censor for A = 0.
#   Note that the LHS is taken as a Surv() object as defined by package survival.
#
# @param modelA1 A formula object. The coxph model for Censor for A = 1.
#   Note that the LHS is taken as a Surv() object as defined by package survival.
#
# @param data A data.frame object. All data required for the analysis.
#
# @param dObj A list object. The column headers of data containing required
#   variables
#
# @param timesB A vector object. The infection times during unblinded phase.
#
# @param timesU A vector object. The infection times during blinded phase.
#
# @param L A numeric object. The study end point.
#
# returns a list containing
#    "censorFit0" the regression object returned for the A = 0 analysis
#    "censorFit1" the regression object returned for the A = 1 analysis
#    "LambdaTU" {nTimesU x 2} Lamdba(u) on timesU for A = 0/1
#    "LambdaTB" {nTimesB x 2} Lamdba(u) on timesB for A = 0/1
#    "LambdaR" {n+2 x 2} Lambda(r) for A = 0/1 with means tacked on
#    "expXbetaG0" {n+2 x 2} exp(xbeta) when G = 0 for A = 0/1
#    "expXbetaG1" {n+2 x 2} exp(xbeta) when G = 1 for A = 0/1
#
#' @importFrom stats reformulate terms coef
#' @import survival
#'
.censorFit <- function(modelA0, 
                       modelA1,  
                       data,
                       dObj,  
                       timesB,  
                       timesU,
                       L) {

  message("\tcensoring")

  n <- nrow(x = data)

  ### create time dependent covariate G(u) = I(R <= u, Gamma = 1)

  # generate random name for id variable and assign 1:n
  while (TRUE) {
    idName = paste(sample(x = LETTERS, size = 10, replace = TRUE), 
                   collapse = "")
    if (idName %in% colnames(x = data)) next
    break
  }
  data[,idName] <- 1L:n

  # prepare base data setting stop time to be U
  args <- list()
  args[[ "data1" ]] <- quote(expr = data)
  args[[ "data2" ]] <- quote(expr = data)
  args[[ "id" ]] <- str2lang(s = idName)
  args[[ "tstop" ]] <- quote(expr = U)

  step1Result <- do.call(what = survival::tmerge, args = args)

  # generate random name for G(u) variable
  while (TRUE) {
    gName = paste(sample(x = LETTERS, size = 10, replace = TRUE), 
                  collapse = "")
    if (gName %in% colnames(x = data)) next
    break
  }

  gam1 <- data[,dObj$Gam] == 1L

  # for each of the Gam==1 subset, covariate "G(u)" is 0 for u in [0,R) and 
  # 1 for u in [R,L)
  # first generate data only for those that experience a transition in
  # the value of this covariate
  nGam1 <- sum(gam1)
  df0 <- data.frame(data[gam1,idName],
                    rep(x = 0.0, times = nGam1),
                    data[gam1,dObj$R],
                    rep(x = 0L, times = nGam1))
  colnames(x = df0) <- c(idName, "start", "stop", gName)

  # now the transition point
  df1 <- data.frame(data[gam1,idName],
                    data[gam1,dObj$R],
                    rep(x = L, times = nGam1),
                    rep(x = 1L, times = nGam1))
  colnames(x = df1) <- c(idName, "start", "stop", gName)

  # and finally the data for those that do not transition
  # for each of this subset, covariate "G" is 0 for u in [0,L)
  nGam0 <- sum(!gam1)
  df <- data.frame(data[!gam1,idName],
                   rep(x = 0.0, times = nGam0),
                   rep(x = L, times = nGam0),
                   rep(x = 0L, times = nGam0))
  colnames(x = df) <- c(idName, "start", "stop", gName)

  # time-dependent data.frame for G(u)
  df <- rbind(df, df1, df0)

  # merge in G() data to base data.frame
  args <- list()
  args[[ "data1" ]] <- quote(expr = step1Result)
  args[[ "data2" ]] <- quote(expr = df)
  args[[ "id" ]] <- str2lang(s = idName)
  args[[ gName ]] <- str2lang(s = paste0("tdc(start, ", gName, ")"))

  step2Result <- do.call(what = survival::tmerge, args = args)

  # generate random name for status variable
  while (TRUE) {
    statusName = paste(sample(x = LETTERS, size = 10, replace = TRUE), 
                       collapse = "")
    if (statusName %in% colnames(x = data)) next
    break
  }

  # add dummy variable for status and assign to 1 if delta == 0
  # this is necessary because macs don't seem to like the use of 1*{} ~ -1
  step2Result[,statusName] <- as.integer(x = step2Result[,dObj$Delta] == 0L)

  # update user formulae to include Surv() LHS and time-dependent covariate G(u)
  survForm <- paste0("Surv(tstart, tstop, ", statusName, ")")

  modelA0 <- stats::reformulate(termlabels = c(attr(x = terms(x = modelA0), 
                                                    which = "term.labels"), 
                                               gName),
                                response = survForm, 
                                intercept = FALSE)

  modelA1 <- stats::reformulate(termlabels = c(attr(x = terms(x = modelA1), 
                                                    which = "term.labels"), 
                                               gName),
                                response = survForm, 
                                intercept = FALSE)

  #  Fit Cox model for U using only placebo (A = 0) participants
  subsetA0 <- step2Result[,dObj$A] == 0L
  A0Fit <- tryCatch(expr = survival::coxph(formula = modelA0, 
                                           data = step2Result,
                                           subset = subsetA0),
                    error = function(e) {
                              stop("unable to fit Cox model for censoring, A = 0\n",
                                   e$message, call. = FALSE)
                            })

  if (anyNA(x = coef(object = A0Fit))) {
    stop("fit of Cox model for censoring, A = 0 resulted in NA coefficients",
         call. = FALSE)
  }

  #  Fit Cox model for U using only placebo (A = 1) participants
  subsetA1 <- step2Result[,dObj$A] == 1L
  A1Fit <- tryCatch(expr = survival::coxph(formula = modelA1, 
                                           data = step2Result,
                                           subset = subsetA1),
                    error = function(e) {
                              stop("unable to fit Cox model for censoring, A = 1\n",
                                   e$message, call. = FALSE)
                            })

  if (anyNA(x = coef(object = A1Fit))) {
    stop("fit of Cox model for censoring, A = 1 resulted in NA coefficients",
         call. = FALSE)
  }

  # Lambda on infection times during blinded and unblinded stages and on R
  # note that lambdaR will have two extra elements tacked on to
  # simplify later Rcpp logic

  # warnings are suppressed because survfit squacks when interactions are used
  # and a data.frame is not provided

  if (!is.null(x = timesB)) {

    LambdaTBA0 <- suppressWarnings(expr = summary(object = survfit(formula = A0Fit),
                                                  time = timesB,
                                                  extend = TRUE)$cumhaz)
    LambdaTBA1 <- suppressWarnings(expr = summary(object = survfit(formula = A1Fit),
                                                  time = timesB,
                                                  extend = TRUE)$cumhaz)
  } else {
    LambdaTBA0 <- 0.0
    LambdaTBA1 <- 0.0
  }

  if (!is.null(x = timesU)) {
    LambdaTUA0 <- suppressWarnings(expr = summary(object = survfit(formula = A0Fit),
                                                  time = timesU,
                                                  extend = TRUE)$cumhaz)
    LambdaTUA1 <- suppressWarnings(expr = summary(object = survfit(formula = A1Fit),
                                                  time = timesU,
                                                  extend = TRUE)$cumhaz)
  } else {
    LambdaTUA0 <- 0.0
    LambdaTUA1 <- 0.0
  }

  # zeros are tacked on to keep Lambdas of same dimension as expXbeta, which
  # includes two rows for means
  LambdaRA0 <- suppressWarnings(expr = summary(object = survfit(formula = A0Fit),
                                               time = c(data[,dObj$R],0.0,0.0),
                                               extend = TRUE)$cumhaz)
  LambdaRA1 <- suppressWarnings(expr = summary(object = survfit(formula = A1Fit),
                                               time = c(data[,dObj$R],0.0,0.0),
                                               extend = TRUE)$cumhaz)

  # expXbeta for all combinations of A = 0, 1 and G = 0, 1
  # exp(Xbeta) includes the average values for each treatment group
  newdata <- data
  newdata <- rbind(newdata,  
                   colMeans(data[data[,dObj$A] == 0L,,drop=FALSE]), 
                   colMeans(data[data[,dObj$A] == 1L,,drop=FALSE]))

  newdata[,gName] <- 0L
  expXbetaA0G0 <- predict(object = A0Fit, newdata = newdata, type = "risk")
  expXbetaA1G0 <- predict(object = A1Fit, newdata = newdata, type = "risk")

  newdata[,gName] <- 1L
  expXbetaA0G1 <- predict(object = A0Fit, newdata = newdata, type = "risk")
  expXbetaA1G1 <- predict(object = A1Fit, newdata = newdata, type = "risk")

  return( list("censorFit0" = A0Fit,
               "censorFit1" = A1Fit,
               "LambdaTU" = cbind(LambdaTUA0, LambdaTUA1),
               "LambdaTB" = cbind(LambdaTBA0, LambdaTBA1),
               "LambdaR" = cbind(LambdaRA0, LambdaRA1),
               "expXbetaG0" = cbind(expXbetaA0G0, expXbetaA1G0),
               "expXbetaG1" = cbind(expXbetaA0G1, expXbetaA1G1)) )
}
