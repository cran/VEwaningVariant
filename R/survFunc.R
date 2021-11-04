#' @import survival
.survFunc <- function(object, newdata) {

  return( exp(x = -predict(object = object,
                           newdata = newdata,
                           type = "expected")) )
}
