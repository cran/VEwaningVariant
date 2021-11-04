#' @import survival
.pzero <- function(object, newdata) {
  temp <- predict(object = object, 
                  newdata = newdata, 
                  type = "lp", 
                  reference = "sample") + 
          sum(coef(object = object) * object$means, na.rm=TRUE)
  return( exp(x = unname(obj=temp)) )
}

