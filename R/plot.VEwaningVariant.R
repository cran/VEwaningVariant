#' Plot Analysis Results
#'
#' Plot the Estimated Vaccine Efficacy
#'
#' When the variant under analysis is present only in the unblinded phase or
#'  only the unblinded phase was included in the analysis, vaccine efficacy
#'  cannot be estimated. In this case, plot() shows the infection rate at 
#'  \eqn{\tau}{tau} = \eqn{\ell}{l} divided by the infection rate at 
#'  \eqn{\tau}{tau}. Recall that \eqn{\tau}{tau} is the time since vaccination,
#'  and \eqn{\ell}{l} is the period of time required to reach full efficacy
#'  after the initial vaccination dose.
#'
#' @param x An object of class VEwaningVariant. The object returned by a call to
#'   veWaningVariant()
#'
#' @param y Ignored
#'
#' @param ... Ignored
#'
#' @param taus A numeric vector object or NULL. The tau values at which
#'   estimates are to be calculated. If NULL, a vector of length nTau
#'   spanning the range [lag, maxTau], where maxTau is the maximum tau
#'   identified from the training data; if xlim is specified, 
#'   [xlim[1], xlim[2]] is used instead. Note that all values must lie in the
#'   range [lag, maxTau].
#'
#' @param xlim A numeric vector object of length 2 or NULL. The extrema
#'   of the tau values at which estimates
#'   are to be calculated. A vector of length nTau spanning the
#'   range [xlim[1], xlim[2]] is generated. Note that the specified limits must
#'   lie in the range [lag, maxTau]. If taus is specified, this input is ignored.
#'
#' @param nTau An integer object. The number of tau values at which
#'   estimates are provided. The default is 20. If taus
#'   is specified, this input is ignored.
#'
#' @param ylim A numeric vector object or NULL. The y-axis limits for the plots.
#'   If NULL, the y-axis limits are taken from the estimated values.
#'
#' @returns A gg object.
#'
#' @name plot
#' @examples
#' data(variantData)
#'
#' set.seed(1234)
#'
#' ind <- sample(1:nrow(x = variantData), 2500, FALSE)
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
                                 taus = NULL, 
                                 xlim = NULL, 
                                 nTau = 20L, 
                                 ylim = NULL) {

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

  theta <- x$theta
  cov <- x$cov

  if (type == 1L) {
    # if only unblinded phase included in original analysis, need to pad
    # theta and cov with intercept term
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

    minTau <- taus[1L]
    maxTau <- utils::tail(x = taus, n = 1L)

  } else if (!is.null(x = xlim)) {

    if (length(x = xlim) != 2L) stop("length(xlim) != 2", call. = FALSE)
    if (!is.numeric(x = xlim)) stop("xlim must be numeric", call. = FALSE)

    # if provided limits, create taus
    xlim <- sort(x = xlim)
    minTau <- max(xlim[1L], lag)
    maxTau <- min(xlim[2L], maxTau)
    taus <- seq(from = minTau, to = maxTau, length.out = nTau)
  } else {
    minTau <- lag
    maxTau <- maxTau
    taus <- seq(from = minTau, to = maxTau, length.out = nTau)
  }

  # limits for generated plot
  xlim <- c(minTau, maxTau)

  # times at which g(u; theta) is calculated are u = tau - lag
  times <- taus - lag

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
    title <- expression(paste("VE(", tau, ")"))
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

  df <- data.frame(taus, veLow, ve, veHigh)

  ggplot(df, aes(x = taus, y = ve)) +
    geom_path(lwd = 1) +
    geom_ribbon(aes(ymin = veLow, ymax = veHigh), alpha = 0.2, fill = "grey") +
    labs(y = title, x = expression(tau)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    coord_cartesian(xlim = xlim, ylim = ylim)

}
