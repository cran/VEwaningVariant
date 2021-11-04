# @param data A data.frame object. All data required for the analysis.
#
# @param model A formula object. The logistic model for Psi.
#
# @param dObj A list object. The column headers of data containing required
#   variables
#
#  returns a vector {n} p_Psi
#
#' @importFrom stats terms reformulate glm predict.glm coef
.psiFit <- function(data, model, dObj) {

  message("\tPsi")

  # use only data for placebo participants (A = 0) with Gamma = 1
  # (have been unblinded)
  subset <- {data[,dObj$A] == 0L} & {data[,dObj$Gam] == 1L}
  means <- colMeans(x = data[subset,])

  # update formula to include Psi
  model <- stats::reformulate(termlabels = attr(x = terms(x = model), 
                                                which = "term.labels"),
                              response = dObj$Psi)

  # fit logistic regression model using only subset of data
  logist <- tryCatch(expr = stats::glm(formula = model,
                                       data = data[subset,],
                                       family = 'binomial'),
                     error = function(e) {
                               stop("unable to obtain fit for Psi\n",
                                    e$message, call. = FALSE)
                             })
  if (anyNA(x = coef(object = logist))) {
    stop("glm fit for Psi returns NA coefficients", call. = FALSE)
  }

  ##  Get predicted probabilities for Psi = 1 for all participants  
  pPsi <- stats::predict.glm(object = logist, newdata = data, type = "response")

  ##  Get predicted probabilities at the means for Psi = 1
  pPsi.mean <- stats::predict.glm(object = logist, 
                                  newdata = as.data.frame(x = t(x = means)), 
                                  type = "response")

  ##  Get stabilized probabilities
    
  pPsi.stab <- pPsi.mean / pPsi

  return( pPsi.stab )
   
}
