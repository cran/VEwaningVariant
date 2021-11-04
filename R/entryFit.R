##  Fit a Cox model to the entry times as a function of X and return a
##  n x 1 vector of stabilized densities evaluated at the observed E, X

# @param data A data.frame object. All data required for the analysis.
#
# @param model A formula object. The coxph model for Entry time. Note
#   that the LHS is taken to be a Surv() object as defined by package survival.
#
# @param dObj A list object. The column headers of data containing required
#   variables
#
#  returns a vector {n} f_E
#
#' @importFrom stats reformulate terms coef
#' @import survival
#
#' @include survFunc.R pzero.R
.entryFit <- function(data, model, dObj) {

  message("\tentry time")

  # create internal naming convention for Status
  while (TRUE) {
    statusName = paste(sample(x = letters, size = 10, replace = TRUE), 
                       collapse = "")
    if (statusName %in% colnames(x = data)) next
    break
  }

  # add dummy variable for status -- everyone in data clearly entered the study
  data[,statusName] <- 1L

  # update formula with Surv() object in LHS
  entryForm <- paste0('Surv(', dObj$E, ',', statusName,')')
  model <- stats::reformulate(termlabels = attr(x = terms(x = model), 
                                                which = "term.labels"),
                              response = entryForm, 
                              intercept = FALSE)

  #  Fit Cox model for E using all participants
  fitObj <- tryCatch(expr = survival::coxph(formula = model, data = data),
                     error = function(e) {
                               stop("unable to fit entry time Cox model\n",
                                    e$message, call. = FALSE)
                             })

  if (anyNA(x = coef(object = fitObj))) {
    stop("fit of Cox model for entry time resulted in NA coefficients",
         call. = FALSE)
  }

  ##  Predicted survival probabilities at each E

  # {n}
  SE <- .survFunc(object = fitObj, newdata = data)

  # {n}
  eE <- .pzero(object = fitObj, newdata = data)

  ##  Predicted survival probabilities at each E for mean of X    

  newdata <- data
  means <- colMeans(x = data)
  for (i in 1L:ncol(x = newdata)) {
    newdata[,i] <- means[i]
  }
  newdata[,dObj$E] <- data[,dObj$E]
  newdata[,statusName] <- 1L
  
  # {n}
  SE.mean <- .survFunc(object = fitObj, newdata = newdata)

  # {n}
  eE.mean <- .pzero(object = fitObj, newdata = newdata)

  ##  Stabilized density for each individual

  # {n}
  fE.stab <- {eE.mean*SE.mean} / {eE*SE}
    
  return( fE.stab )
}
