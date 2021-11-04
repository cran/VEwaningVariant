setGeneric(name = ".verifyModels",
           def = function(entry, 
                          unblind0, 
                          unblind1, 
                          psi, 
                          censor0, 
                          censor1, ...) { 
               standardGeneric(".verifyModels") 
             })

setMethod(f = ".verifyModels",
          signature = c(entry = "ANY", 
                        unblind0 = "ANY", 
                        unblind1 = "ANY", 
                        psi = "ANY", 
                        censor0 = "ANY", 
                        censor1 = "ANY"),
          definition = function(entry, 
                                unblind0, 
                                unblind1, 
                                psi, 
                                censor0, 
                                censor1, ...) { 
              stop("inappropriate combination of models specified for weights", 
                   call. = FALSE)
            })

setMethod(f = ".verifyModels",
          signature = c(entry = "NULL", 
                        unblind0 = "NULL", 
                        unblind1 = "NULL", 
                        psi = "NULL", 
                        censor0 = "NULL", 
                        censor1 = "NULL"),
          definition = function(entry, 
                                unblind0, 
                                unblind1, 
                                psi, 
                                censor0, 
                                censor1, ...) { 
              message("no models provided for weights -- weights = 1")
              return( list("covs" = NULL, "wgtType" = 0L) )
            })

setMethod(f = ".verifyModels",
          signature = c(entry = "NULL", 
                        unblind0 = "NULL", 
                        unblind1 = "NULL", 
                        psi = "NULL", 
                        censor0 = "formula", 
                        censor1 = "formula"),
          definition = function(entry, 
                                unblind0, 
                                unblind1, 
                                psi, 
                                censor0, 
                                censor1, ...) { 

              message("weighting depends only on censoring dynamic")

              covs <- .extractCov(model = censor0)
              covs <- c(covs, .extractCov(model = censor1))

              # keep only 1 copy of covariates
              covs <- sort(x = unique(x = covs))

              return( list("covs" = covs, "wgtType" = 1L) )
            })


setMethod(f = ".verifyModels",
          signature = c(entry = "formula", 
                        unblind0 = "formula", 
                        unblind1 = "formula", 
                        psi = "formula", 
                        censor0 = "formula", 
                        censor1 = "formula"),
          definition = function(entry, 
                                unblind0, 
                                unblind1, 
                                psi, 
                                censor0, 
                                censor1, ...) { 

              message("models provided for all components of weights")

              covs <- .extractCov(model = entry)
              covs <- c(covs, .extractCov(model = unblind0))
              covs <- c(covs, .extractCov(model = unblind1))
              covs <- c(covs, .extractCov(model = psi))
              covs <- c(covs, .extractCov(model = censor0))
              covs <- c(covs, .extractCov(model = censor1))

              # keep only 1 copy of covariates
              covs <- sort(x = unique(x = covs))

              return( list("covs" = covs, "wgtType" = 2L) )
            })

.extractCov <- function(model) {

  # remove possible LHS of provided model
  model <- stats::update.formula(old = model, new = NULL ~ .)

  # extract covariate names from factors attribute of terms object
  covs <- rownames(x = attr(x = stats::terms(x = model), which = "factors"))

  return( covs )
}
