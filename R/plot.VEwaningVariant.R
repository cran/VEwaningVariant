#' Plot Analysis Results
#'
#' Plot the Estimated Vaccine Efficacy
#'
#' When the variant under analysis is present only in the unblinded phase, 
#'   vaccine efficacy cannot be estimated. In this case, plot() shows the  
#'   relative infection rate at times t since full efficacy reached, defined  
#'   as infection rate at time t = time since full efficacy reached  
#'   divided by the infection rate at the time full efficacy is reached (t=0).
#'
#' @param x An object of class VEwaningVariant. The object returned by a call to
#'   veWaningVariant()
#'
#' @param y Ignored
#'
#' @param ... Ignored
#'
#' @param times A numeric vector object or NULL. The times since full
#'   efficacy at which the vaccine efficacy is to be estimated. If NULL, the
#'   times will be generated internally as a vector of length nTimes spanning
#'   the range [0, maxTime], where maxTime is the maximum time since vaccination  
#'   present in the original analysis. Values provided outside of [0, maxTime] 
#'   are ignored.
#'
#' @param xlim A numeric vector object of length 2 or NULL. The extrema
#'   of the times values at which estimates are to be calculated. A vector of 
#'   length nTimes spanning the range [xlim[1], xlim[2]] is generated. Note 
#'   that the specified limits must lie in the range [0, maxTime]. If input 
#'   times is a vector object, this input is ignored.
#'
#' @param nTimes An integer object. The number of time values at which
#'   estimates are obatined. The default is 20. If input times is a vector
#'   object, this input is ignored.
#'
#' @param ylim A numeric vector object or NULL. The y-axis limits for the plots.
#'   If NULL, the y-axis limits are taken from the estimated values.
#'
#' @returns A gg object.
#'
#' @name plot
#' @examples
#' \dontshow{
#'   RcppArmadillo::armadillo_throttle_cores(2)
#' }
#' data(variantData)
#'
#' set.seed(1234)
#'
#' ind <- sample(1:nrow(x = variantData), 2000, FALSE)
#'
#' # NOTE: This sample size is chosen for example only -- larger data sets
#' # should be used.
#'
#' res <- veWaningVariant(data = variantData[ind,], 
#'                        L = 52,  
#'                        lag = 6,
#'                        gFunc = 'piece', 
#'                        v = c(15,30))
#'
#' plot(x = res)
#' @method plot VEwaningVariant
#' @export 
#' @importFrom graphics plot
#' @import ggplot2
#' @importFrom utils tail
#
plot.VEwaningVariant <- function(x, y, ..., 
                                 times = NULL, 
                                 xlim = NULL, 
                                 nTimes = 20L, 
                                 ylim = NULL) {

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

  theta <- x$theta
  cov <- x$cov

  if (type == 1L) {
    # if only unblinded phase included in original analysis, need to pad
    # theta and cov with intercept term
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

    minTime <- times[1L]
    maxTime <- utils::tail(x = times, n = 1L)

  } else if (!is.null(x = xlim)) {

    if (length(x = xlim) != 2L) stop("length(xlim) != 2", call. = FALSE)
    if (!is.numeric(x = xlim)) stop("xlim must be numeric", call. = FALSE)

    # if provided limits, create times
    xlim <- sort(x = xlim)
    minTime <- max(xlim[1L], 0.0)
    maxTime <- min(xlim[2L], maxTime)
    times <- seq(from = minTime, to = maxTime, length.out = nTimes)
  } else {
    minTime <- 0.0
    times <- seq(from = minTime, to = maxTime, length.out = nTimes)
  }

  # limits for generated plot
  xlim <- c(minTime, maxTime)

  gFuncR <- gFunction(gFunc = gFunc, u = times, theta = theta, knots = v)

  dg <- matrix(data = gFuncR[[ 2L ]], ncol = length(x = theta))

  se <- sqrt(x = diag(x = dg %*% cov %*% t(x = dg)))

  if (type == 1L) {
    # theta0 is not estimated when only unblinded phase is included
    veHigh <- exp(x = -{gFuncR[[ 1L ]] - 1.96*se})
    ve <- exp(x = -gFuncR[[ 1L ]])
    veLow <- exp(x = -{gFuncR[[ 1L ]] + 1.96*se})
    title <- "Relative Infection Rate"
  } else {
    veLow <- 1.0 - exp(x = theta[1L] + gFuncR[[ 1L ]] - 1.96*se)
    ve <- 1.0 - exp(x = theta[1L] + gFuncR[[ 1L ]])
    veHigh <- 1.0 - exp(x = theta[1L] + gFuncR[[ 1L ]] + 1.96*se)
    title <- "VE(t)"
  }

  # limits of the y-axis for generated plots
  if (is.null(x = ylim)) {
    ylim <- round(c(min(veLow, ve, veHigh) - 0.1, 
                    max(veLow, ve, veHigh) + 0.1),2L)
  } else {
    if (length(x = ylim) != 2L) stop("length(ylim) != 2", call. = FALSE)
    if (!is.numeric(x = ylim)) stop("ylim must be numeric", call. = FALSE)
    ylim <- sort(x = ylim)
  }

  df <- data.frame(times, veLow, ve, veHigh)

  ggplot(df, aes(x = times, y = ve)) +
    geom_path(lwd = 1) +
    geom_ribbon(aes(ymin = veLow, ymax = veHigh), alpha = 0.2, fill = "grey") +
    labs(y = title, x = "Time since full efficacy (t)") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    coord_cartesian(xlim = xlim, ylim = ylim)

}
