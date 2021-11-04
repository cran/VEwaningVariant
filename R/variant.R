# identify participants to include in the analysis based on the variant
# value provided by the user

setGeneric(name = ".variant",
           def = function(variant, ...) { 
               standardGeneric(".variant") 
             })

# anything not explicitly allowed is forbidden
setMethod(f = ".variant",
          signature = c(variant = "ANY"),
          definition = function(variant, ...) { 
              stop("inappropriate object provided for variant", call. = FALSE)
            })

# NULL is taken as the "all variant" case -- though this case is expected
# to be provided as 0 by the user
setMethod(f = ".variant",
          signature = c(variant = "NULL"),
          definition = function(variant, delta, ...) { 

              message("all variant infections included in analysis")

              group <- rep(x = 1L, times = length(x = delta))

              # remove cases that did not experience infection
              group[delta <= 0L] <- 0L

              return( list("group" = group) )
            })

setMethod(f = ".variant",
          signature = c(variant = "integer"),
          definition = function(variant, delta, ...) { 

              if (variant > 0L) {
                # variant specific analysis
                group <- rep(x = 0L, times = length(x = delta))

                group[delta == variant] <- 1L

                message("only variant infections ", variant, 
                        " included in analysis")

                ng <- sum(group)
                ni <- sum(delta > 0L)

                if (ng == 0L) {
                  stop("no variant ", variant, " infections found in data",
                       call. = FALSE)
                } else if (ng < ni*0.1) {
                  message("NOTE: < 10% of infections are variant ", variant)
                }

                return( list("group" = group) )

              } else if (variant == 0L) {

                # if all variants, call NULL method
                return( .variant(variant = NULL, delta = delta) )

              } else {
                stop("inappropriate value provided for variant", call. = FALSE)
              }

            })

# if passed as a non-integer numeric, recast as integer and call integer method 
setMethod(f = ".variant",
          signature = c(variant = "numeric"),
          definition = function(variant, ...) { 
              iv <- as.integer(round(x = variant, digits = 0L))
              if (!isTRUE(all.equal(iv, variant))) {
                stop("variant must be integer", call. = FALSE)
              }
              return( .variant(variant = iv, ...) )
            })
