#' @useDynLib VEwaningVariant, .registration=TRUE
#' @include RcppExports.R

esttheta <- function(data, 
                     timesB, 
                     timesU,
                     fitERPsi, 
                     censorFit, 
                     gFunc,
                     dObj,
                     minWgt,
                     maxWgt,
                     wgtType,
                     phaseType) {

  message("estimating theta")

  # Rcpp doesn't like NULLs
  if (is.null(x = timesB)) timesB <- 0.0
  if (is.null(x = timesU)) timesU <- 0.0

  ##  Initialize -- starting value 

  theta <- rep(x = 0.0, times = gFunc$nTheta)
  score <- rep(x = 1.0, times = gFunc$nTheta)

  bmax <- Inf
  tol <- 1e-5
  imax <- 20L
  iter <- 1L

  message("iteration:")  

  while (iter < imax && bmax > tol) {

    message(sprintf("%3d",iter), appendLF = {iter == 10L})  

    eu <- estTheta(E = data[,dObj$E],
                   U = data[,dObj$U],
                   R = data[,dObj$R],
                   lag = data[,dObj$lag],
                   A = data[,dObj$A],
                   Gam = data[,dObj$Gam],
                   Psi = data[,dObj$Psi],
                   delta = data[,dObj$Delta],
                   wgt = wgtType,
                   fR = fitERPsi$unblindFit$fR,
                   fEX = fitERPsi$fEX,
                   Psiprob = fitERPsi$Psiprob,
                   minWgt = minWgt,
                   maxWgt = maxWgt,
                   censor_expXbetaG0 = censorFit$expXbetaG0,
                   censor_expXbetaG1 = censorFit$expXbetaG1,
                   censor_LambdaTB = censorFit$LambdaTB,
                   censor_LambdaTU = censorFit$LambdaTU,
                   censor_LambdaR = censorFit$LambdaR,
                   unBlind_expXbeta = fitERPsi$unblindFit$expXbeta,
                   unBlind_LambdaT =  fitERPsi$unblindFit$LambdaT,
                   timesB = timesB,  
                   timesU = timesU,  
                   gFunc = gFunc$gFunc,
                   group = gFunc$group,
                   v = gFunc$v,
                   theta = theta,
                   type = phaseType)

    score <- eu[[ 1L ]]
    grad <- matrix(eu[[ 2L ]], nrow = gFunc$nTheta)

    if (anyNA(x = score) || anyNA(x = grad)) {
      stop("\nencountered NaNs in Newton-Raphson", call. = FALSE)
    }

    ##  Update theta
    if (phaseType == 1L) {
      # if only unblinded phase, cannot estimate theta0
      score <- score[-1L]
      grad <- grad[-1L,,drop=FALSE]
      grad <- grad[,-1L,drop=FALSE]
      inv <- tryCatch(expr = solve(a = grad, b = score),
                      error = function(e) {
                                stop("unable to invert gradiant in Newton-Raphson\n",
                                     e$message, call. = FALSE)
                              })
      theta <- theta[-1L] + inv
      theta <- c(0.0, theta)
    } else {
      inv <- tryCatch(expr = solve(a = grad, b = score),
                      error = function(e) {
                                stop("unable to invert gradiant in Newton-Raphson\n",
                                     e$message, call. = FALSE)
                              })
      theta <- theta + inv
    }

    iter <- iter + 1L
    bmax <- max(abs(x = score))
  }
  message("")

  ##  For final value of theta, get covariance matrix and compute estimated 
  ##  SEs -- we only use sandwich

  eu <- estTheta(E = data[,dObj$E],
                 U = data[,dObj$U],
                 R = data[,dObj$R],
                 lag = data[,dObj$lag],
                 A = data[,dObj$A],
                 Gam = data[,dObj$Gam],
                 Psi = data[,dObj$Psi],
                 delta = data[,dObj$Delta],
                 wgt = wgtType,
                 fR = fitERPsi$unblindFit$fR,
                 fEX = fitERPsi$fEX,
                 Psiprob = fitERPsi$Psiprob,
                 minWgt = minWgt,
                 maxWgt = maxWgt,
                 censor_expXbetaG0 = censorFit$expXbetaG0,
                 censor_expXbetaG1 = censorFit$expXbetaG1,
                 censor_LambdaTB = censorFit$LambdaTB,
                 censor_LambdaTU = censorFit$LambdaTU,
                 censor_LambdaR = censorFit$LambdaR,
                 unBlind_expXbeta = fitERPsi$unblindFit$expXbeta,
                 unBlind_LambdaT =  fitERPsi$unblindFit$LambdaT,
                 timesB = timesB,  
                 timesU = timesU,  
                 gFunc = gFunc$gFunc,
                 group = gFunc$group,
                 v = gFunc$v,
                 theta = theta,
                 type = phaseType)

  score <- eu[[ 1L ]]
  grad <- matrix(eu[[ 2L ]], nrow = gFunc$nTheta)
  meat <- matrix(eu[[ 3L ]], nrow = gFunc$nTheta)

  if (phaseType == 1L) {
    # if only unblinded phase, cannot estimate theta0
    theta <- theta[-1L]
    grad <- grad[-1L,,drop=FALSE]
    grad <- grad[,-1L,drop=FALSE]
    meat <- meat[-1L,,drop=FALSE]
    meat <- meat[,-1L,drop=FALSE]
    Covmodel <- tryCatch(expr = solve(a = grad),
                         error = function(e) {
                                   stop("unable to invert gradiant for sandwich estimator\n",
                                        e$message, call. = FALSE)
                                 })
    nms <- paste0('theta', 1L:{gFunc$nTheta- 1L})
  } else {
    Covmodel <- tryCatch(expr = solve(a = grad),
                         error = function(e) {
                                   stop("unable to invert gradiant for sandwich estimator\n",
                                        e$message, call. = FALSE)
                                 })
    nms <- paste0('theta', 0L:{gFunc$nTheta- 1L})
  }

  Covsand <- Covmodel %*% meat %*% Covmodel
  SEsand <- sqrt(x = diag(x = Covsand))

  names(x = theta) <- nms

  colnames(x = Covsand) <- nms
  rownames(x = Covsand) <- nms
  names(x = SEsand) <- nms

  return( list("theta" = theta, "cov" = Covsand, "SE" = SEsand) )

}
