# verify appropriate specification of gFunc and process nodes
#
# piece = 1; linear = 2; splin = 3; spcub = 4

setGeneric(name = ".gFunc_v",
           def = function(gFunc, v, ...) { 
               standardGeneric(".gFunc_v") 
             })

# anything not explicitly allowed is forbidden
setMethod(f = ".gFunc_v",
          signature = c(gFunc = "ANY",
                        v = "ANY"),
          definition = function(gFunc, v, ...) { 
              stop("inappropriate object(s) provided for gFunc and/or v", 
                   call. = FALSE)
            })

# If gFunc is NULL, g(theta) is taken as a linear function
setMethod(f = ".gFunc_v",
          signature = c(gFunc = "NULL",
                        v = "ANY"),
          definition = function(gFunc, v, ...) { 
              # default g-function is linear, v is set to 1 for Rcpp
              return( list("gFunc" = 2L, "v" = 1.0, "nTheta" = 2L) )
            })

# if gFunc is a character and no nodes are provided, gFunc must be 'lin'
setMethod(f = ".gFunc_v",
          signature = c(gFunc = "character",
                        v = "NULL"),
          definition = function(gFunc, v, ...) { 

              if (gFunc %in% c("splin", "spcub", "piece")) {
                stop("v must be speficied for gFunc = ", gFunc, call. = FALSE)
              } else if (gFunc != 'lin') {
                stop("gFunc not recognized", call. = FALSE)
              }

              return( .gFunc_v(gFunc = NULL, v = v) )

            })

# If character, it must be 'piece', 'lin', 'splin', or 'spcub'
setMethod(f = ".gFunc_v",
          signature = c(gFunc = "character",
                        v = "numeric"),
          definition = function(gFunc, v, L, ...) { 

              # ensure that knots are unique and sorted
              v <- sort(x = unique(x = v))

              if (gFunc == "lin") {
		# if 'lin', call NULL method
                return( .gFunc_v(gFunc = NULL, v = v) )

              } else if (gFunc == "piece") {

                # piecewise needs lower and upper bounds to be 0, L respectively
                v <- unique(x = c(0.0, v, L))

                # use only nodes that lie in [0,L]
                v <- v[{v <= {L+1e-8}} & {v > {0.0-1e-8}}]

                if (length(x = v) <= 2L) {
                  stop("inappropriate value provided for v", call. = FALSE)
                }

                return( list("gFunc" = 1L, 
                             "v" = v, 
                             "nTheta" = length(x = v) - 1L) )

              } else if (gFunc == "splin") {

                # use only nodes that lie in (0,L)
                v <- v[{v < L} & {v > 0.0}]

                if (length(x = v) == 0L) {
                  stop("inappropriate value provided for v", call. = FALSE)
                }

                # linear spline has two additional thetas
                return( list("gFunc" = 3L, 
                             "v" = v, 
                             "nTheta" = length(x = v) + 2L) )

              } else if (gFunc == "spcub") {

                # use only nodes that lie in (0,L)
                v <- v[{v < L} & {v > 0.0}]

                if (length(x = v) == 0L) {
                  stop("inappropriate value provided for v", call. = FALSE)
                }

                # cubic spline has four additional thetas
                return( list("gFunc" = 4L, 
                             "v" = v, 
                             "nTheta" = length(x = v) + 4L) )

              } else {
                stop("gFunc not recognized", call. = FALSE)
              }


            })
