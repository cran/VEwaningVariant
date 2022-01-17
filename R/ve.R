#' Retrieve the Estimated Vaccine Efficacy
#'
#' Uses a prior veWaningVariant() analysis to estimate the vaccine efficacy
#'   at the provided times since full efficacy.
#'
#' When the variant under analysis is present only in the unblinded phase, 
#'   vaccine efficacy cannot be estimated. In this case, ve() returns the  
#'   relative infection rate at times t since full efficacy reached, defined  
#'   as infection rate at time t = time since full efficacy reached  
#'   divided by the infection rate at the time full efficacy is reached (t=0).
#'
#' @param x An object of class VEwaningVariant. The object returned by a call to
#'   veWaningVariant()
#'
#' @param times A numeric vector object or NULL. The times since full
#'   efficacy at which the vaccine efficacy is to be estimated. If NULL, the
#'   times will be generated internally as a vector of length nTimes spanning
#'   the range [0, maxTime], where maxTime is the maximum time since vaccination  
#'   present in the original analysis. Values provided outside of [0, maxTime] 
#'   are ignored.
#'
#' @param nTimes An integer object. The number of time values at which
#'   estimates are obtained. The default is 20. If input times is a vector
#'   object, this input is ignored.
#'
#' @returns A matrix object. The first column contains the times since
#'   full efficacy at which the estimates are provided; the second column
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
#' ind <- sample(1:nrow(variantData), 2000)
#' # NOTE: This sample size is chosen for example only -- larger data sets
#' # should be used.
#' 
#' res <- veWaningVariant(data = variantData[ind,], 
#'                        L = 52,  
#'                        gFunc = 'piece', 
#'                        v = c(15,30))
#'
#' ve(x = res, times = c(10,20,30,40,50))
#' @export 
ve <- function(x, times = NULL, nTimes = 20L) {

  # attributes contained in the VEwaningVariant object and needed for
  # estimation
  maxTime <- attr(x = x, which = "maxTime")
  gFunc <- attr(x = x, which = "gFunc")
  v <- attr(x = x, which = "v")
  if (is.null(x = v)) v <- 1.0
  type <- attr(x = x, which = "phaseType")

  if (is.null(x = maxTime)) {
    stop("unable to use post-processing tools, ",
         "no participant reached full efficacy in original analysis")
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

  if (!is.integer(x = nTimes)) nTimes <- as.integer(x = nTimes)
  if (nTimes <= 0L) stop("inappropriate value provided for nTimes", call. = FALSE)

  if (!is.null(x = times)) {

    if (!is.numeric(x = times)) stop("times must be numeric", call. = FALSE)

    times <- sort(x = unique(x = times))

    if (times[1L] < 0 || any(times > maxTime)) {
      message("times values outside of [0, maxTime] are ignored")
      times <- times[{times > -1e-8} & {times < {maxTime+1e-8}}]
      if (length(x = times) == 0L) {
        stop("inappropriate times values provided", call. = FALSE)
      }
    }

  } else {
    times <- seq(from = 0.0, to = maxTime, length.out = nTimes)
  }

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
    return( cbind("time" = times, "RelInfRate" = rate, "SE" = se) )
  } else {
    return( cbind("time" = times, "VE" = 1.0 - rate, "SE" = se) )
  }

}
