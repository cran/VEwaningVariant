#' Print Analysis Results
#'
#' Print the primary results of the analysis
#'
#' @param x An object of class VEwaningVariant. The object returned by a call to
#'   veWaningVariant()
#'
#' @param ... Ignored
#'
#'
#' @name print
#'
#' @returns No return value, called to display key results.
#'
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
#' print(x = res)
#' @method print VEwaningVariant
#' @export 
print.VEwaningVariant <- function(x, ...) {

  attr(x = x, which = "gFunc") <- NULL
  attr(x = x, which = "maxTau") <- NULL
  attr(x = x, which = "lag") <- NULL
  attr(x = x, which = "v") <- NULL
  attr(x = x, which = "phaseType") <- NULL
  attr(x = x, which = "wgtType") <- NULL

  x <- unclass(x = x)

  print(x)

  return( NULL )

}
