#' @include confirmIntegerLike.R variant.R gFunc_v.R
#' @importFrom stats complete.cases
.verifyInputs <- function(data, dObj, covs, gFunc, variant, v, L) {

  message("verifying inputs")

  # data must be a data.frame or matrix with column headers
  if (is.matrix(x = data)) {
    if (is.null(x = colnames(x = data))) {
      stop("data must be a data.frame object", call. = FALSE)
    }
    data <- as.data.frame(x = data)
  }
  if (!is.data.frame(x = data)) {
      stop("data must be a data.frame object", call. = FALSE)
  }

  # ensure that the names exist in the provided data.frame

  tryCatch(expr = data[,c(dObj$Gam, dObj$R, dObj$E, dObj$U, 
                          dObj$A, dObj$Psi, dObj$Delta, dObj$lag)],
           error = function(e) {
                     stop("unable to identify key components of input data\n",
                          e$message, call. = FALSE)
                   })

  if (!is.null(x = covs)) {
    tryCatch(expr = data[,covs],
             error = function(e) {
                       stop("unable to identify model covariates in input data\n",
                            e$message, call. = FALSE)
                     })
  }

  # limit data to only those covariates required for analysis
  # this avoids removal of cases with NA in unused data
  data <- data[,c(unlist(x = dObj), covs),drop=FALSE]

  # data must be complete
  complete <- stats::complete.cases(data)
  nRm <- sum(!complete)
  if (nRm > 0L) {
    data <- data[complete,,drop=FALSE]
    message("\tremoved ", nRm, " cases due to incomplete data")
  }

  ## Entry Times

  # ensure that entry times are non-negative
  if (any(data[,dObj$E] < 0.0)) {
    stop("entry time must be non-negative", call. = FALSE)
  }

  ## A

  # ensure that treatment is an integer
  iA <- .confirmIntegerLike(x = data[,dObj$A], name = "treatment")

  # ensure that treatment is one of (0,1)
  if (any(!{iA %in% c(0L,1L)})) {
    stop("unrecognized treatment values", call. = FALSE)
  } else if (all(iA == 0L) || all(iA == 1L)) {
    stop("all participants received the same treatment", call. = FALSE)
  }

  data[,dObj$A] <- iA

  ## U

  # reset U > L and Delta to 0
  tst <- data[,dObj$U] > L
  if (any(tst)) {
    message("\t", sum(tst), 
            " records with U > L; reset cases as U = L and Delta = 0")
  }
  data[tst,dObj$U] <- L
  data[tst,dObj$Delta] <- 0L

  # U must be >= E
  tst <- data[,dObj$U] < data[,dObj$E]
  if (any(tst)) {
    message("\tremoved ", sum(tst), " records with U < E")
    data <- data[!tst,,drop=FALSE]
  }


  ## R

  # If R > L and Gamma = 1 -- unblinding happened first and occurred
  #  after study end point so participant was censored
  # If R > L and Gamma = 0 -- R was the time of infection or censoring
  #   if infection, we've already set U = L and reset Gamma to indicate
  #   no infection, so L is time of censoring; if censoring, just shifted
  #   censoring time to time of study end point
  tst <- data[,dObj$R] > L
  if (any(tst)) {
    message("\t", sum(tst), 
            " records with R > L; reset cases as R = L and Gamma = 0")
  }
  data[tst,dObj$R] <- L
  data[tst,dObj$Gam] <- 0L

  # R must be >= E
  tst <- data[,dObj$R] < data[,dObj$E]
  if (any(tst)) {
    message("\tremoved ", sum(tst), " records with R < E")
    data <- data[!tst,,drop=FALSE]
  }

  # R must be <= U
  tst <- data[,dObj$R] > data[,dObj$U]
  if (any(tst)) {
    message("\tremoved ", sum(tst), " records with R > U")
    data <- data[!tst,,drop=FALSE]
  }

  ## Gamma

  # ensure that the Gamma is an integer
  iGam <- .confirmIntegerLike(x = data[,dObj$Gam], name = "Gamma")

  # ensure that Gamma is one of (0,1)
  if (any(!{iGam %in% c(0L,1L)})) {
    stop("unrecognized Gamma values", call. = FALSE)
  }
  data[,dObj$Gam] <- iGam

  ## Psi

  # ensure that Psi is integer
  iPsi <- .confirmIntegerLike(x = data[,dObj$Psi], name = "Psi")

  # ensure that Psi is one of (0,1)
  if (any(!{iPsi %in% c(0L,1L)})) {
    stop("unrecognized Psi values", call. = FALSE)
  }
  data[,dObj$Psi] <- iPsi

  refused <- {data[,dObj$A] == 0L} & 
             {data[,dObj$Gam] == 1L} & 
             {data[,dObj$Psi] == 0L}
  accepted <- {data[,dObj$A] == 0L} & 
              {data[,dObj$Gam] == 1L} & 
              {data[,dObj$Psi] == 1L}

  if (sum(refused) > 0L & sum(accepted) > 0L) {
    message("\t", sprintf("%10d", sum(refused)), 
            " placebo participants refused vaccine after unblinding\n",
            "\t", sprintf("%10d", sum(accepted)), 
            " placebo participants accepted vaccine after unblinding")
  } else if (sum(refused) == 0L) {
    message("\tall placebo participants accepted vaccine after unblinding")
  } else if (sum(accepted) == 0L) {
    message("\tno placebo participants accepted vaccine after unblinding")
  }

  ## lag
  # again, this is hear only in anticipation of individual specific lag

  # lag cannot be negative
  if (any(data[,dObj$lag] < 0)) {
    stop("lag must be non-negative", call. = FALSE)
  }

  # if provided as infinity, set to value > L
  if (any(is.infinite(x = data[,dObj$lag]))) {
    tst <- is.infinite(x = data[,dObj$lag])
    data[tst,dObj$lag] <- L + 10.0
  }

  ## Delta

  # ensure that variant is an integer
  iv <- .confirmIntegerLike(x = data[,dObj$Delta], name = "Delta")

  # ensure that all are >=-1
  if (any(iv < -1L)) stop("Delta must be >= -1", call. = FALSE)
  
  data[,dObj$Delta] <- iv

  div <- table(iv, useNA = "no")

  for (i in 1L:length(x = div)) {
    if (names(x = div)[i] == "0") {
      message("\t", sprintf("%10d",div[i]), 
              ifelse(test = div[i] > 1, 
                     yes = " participants", 
                     no = " participant"),
              " were censored")
    } else if (names(x = div)[i] == "-1") {
      if (div[i] != sum(refused)) {
        stop("verify data -- the number of participants with ",
             "A = 0, Psi = 0, and Gam = 1 ",
             "does not agree with the number of Delta = -1")
      }
    } else {
      message("\t", sprintf("%10d",div[i]), 
              ifelse(test = div[i] > 1, 
                     yes = " participants", 
                     no = " participant"),
              " experienced variant ", names(x = div)[i], " infection")
    }
  }



  ## gFunc / v

  # ensure gFunc and v are appropriately specified and identify participants
  # that experienced the variant
  gFuncObj <- .variant(variant = variant, delta = data[,dObj$Delta])

  gFuncObj <- c(gFuncObj, .gFunc_v(gFunc = gFunc, v = v, L = L))

  return( list("data" = data, "gFuncObj" = gFuncObj) )

}

