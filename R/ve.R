#' Retrieve the Estimated Vaccine Efficacy
#'
#' Uses a prior veWaningVariant() analysis to estimate the vaccine efficacy
#'   at the provided times since first vaccination dose.
#'
#' When the variant under analysis is present only in the unblinded phase, 
#'  vaccine efficacy
#'  cannot be estimated. In this case, ve() returns the infection rate at 
#'  \eqn{\tau}{tau} = \eqn{\ell}{l} divided by the infection rate at 
#'  \eqn{\tau}{tau}. Recall that \eqn{\tau}{tau} is the time since vaccination,
#'  and \eqn{\ell}{l} is the period of time required to reach full efficacy
#'  after first vaccination dose.
#'
#' @param x An object of class VEwaningVariant. The object returned by a call to
#'   veWaningVariant()
#'
#' @param taus A numeric vector object or NULL. The times since first
#'   vaccination dose at which
#'   the vaccine efficacy is to be estimated. If NULL, a vector of length nTau
#'   spanning the range [lag, maxTau], where maxTau is the maximum tau
#'   identified from the original analysis. Values provided outside of
#'   [lag, maxTau] are ignored.
#'
#' @param nTau An integer object. The number of tau values at which
#'   estimates are provided. The default is 20. If input taus
#'   is specified, this input is ignored.
#'
#' @returns A matrix object. The first column contains the times since
#'   vaccination at which the estimates are provided; the second column
#'   contains estimated vaccine efficacy or relative infection rate 
#'   (see Details); and the third is the standard error.
#'
#'
#' @name ve
#' @examples
#' data(variantData)
#'
#' set.seed(1234)
#'
#' ind <- sample(1:nrow(variantData), 2500)
#' # NOTE: This sample size is chosen for example only -- larger data sets
#' # should be used.
#' 
#' res <- veWaningVariant(data = variantData[ind,], 
#'                        L = 52,  
#'                        gFunc = 'piece', 
#'                        v = c(15,30))
#'
#' ve(x = res, taus = c(10,20,30,40,50))
#' @export 
ve <- function(x, taus = NULL, nTau = 20L) {

  # attributes contained in the VEwaningVariant object and needed for
  # estimation
  lag <- attr(x = x, which = "lag")
  maxTau <- attr(x = x, which = "maxTau")
  gFunc <- attr(x = x, which = "gFunc")
  v <- attr(x = x, which = "v")
  if (is.null(x = v)) v <- 1.0
  type <- attr(x = x, which = "phaseType")

  if (is.null(x = maxTau)) {
    stop("unable to use post-processing tools, all tau < lag in original analysis")
  }

  # if only unblinded phase included in original analysis, need to pad
  # theta and cov with intercept term
  theta <- x$theta
  cov <- x$cov

  # type == 1 indicates that unblinded phase not included in the analysis,
  # must pad theta and cov with theta0
  if (type == 1L) {
    theta <- c(0.0, theta)
    cov <- rbind(0.0, cbind(0.0, cov))
  }

  if (!is.integer(x = nTau)) nTau <- as.integer(x = nTau)
  if (nTau <= 0L) stop("inappropriate value provided for nTau", call. = FALSE)

  if (!is.null(x = taus)) {

    if (!is.numeric(x = taus)) stop("taus must be numeric", call. = FALSE)

    taus <- sort(x = unique(x = taus))

    if (taus[1L] < lag || any(taus > maxTau)) {
      message("tau values outside of [lag, maxTau] are ignored")
      taus <- taus[taus > {lag-1e-8} & taus < {maxTau+1e-8}]
      if (length(x = taus) == 0L) {
        stop("inappropriate tau values provided", call. = FALSE)
      }
    }

  } else {
    taus <- seq(from = lag, to = maxTau, length.out = nTau)
  }

  times <- taus - lag
  gFuncR <- gFunction(gFunc = gFunc, u = times, theta = theta, knots = v)

  if (type == 1L) {
    # type == 1 indicates that blinded phase not included in the analysis,
    # theta0 cannot be estimated
    # return infection rate
    rate <- drop(x = exp(x = -gFuncR[[ 1L ]]))
  } else {
    rate <- drop(x = exp(x = theta[1L] + gFuncR[[ 1L ]]))
  }

  dg <- matrix(data = gFuncR[[ 2L ]], ncol = length(x = theta))

  drate <- dg*rate

  se <- sqrt(x = diag(x = drate %*% cov %*% t(x = drate)))

  if (type == 1L) {
    return( cbind("tau" = taus, "RelInfRate" = rate, "SE" = se) )
  } else {
    return( cbind("tau" = taus, "VE" = 1.0 - rate, "SE" = se) )
  }

}
