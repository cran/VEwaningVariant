setGeneric(name = ".confirmIntegerLike",
           def = function(x, ...) { 
               standardGeneric(".confirmIntegerLike") 
             })

setMethod(f = ".confirmIntegerLike",
          signature = c(x = "ANY"),
          definition = function(x, name, ...) { 
              stop(name, " must be integer-like", call. = FALSE)
            })

# Note this sets the base level to 0
setMethod(f = ".confirmIntegerLike",
          signature = c(x = "factor"),
          definition = function(x, name, ...) { 
              levs <- levels(x = x)
              ix <- match(x = x, table = levs) - 1L
              return( ix )
            })

setMethod(f = ".confirmIntegerLike",
          signature = c(x = "integer"),
          definition = function(x, name, ...) { 
              return( x )
            })

setMethod(f = ".confirmIntegerLike",
          signature = c(x = "numeric"),
          definition = function(x, name, ...) { 
              tx <- as.integer(x = round(x = x, digits = 0L))
              if (!isTRUE(x = all.equal(target = tx, current = x))) {
                stop(name, " must be integer-like", call. = FALSE)
              }
              return( tx )
            })
