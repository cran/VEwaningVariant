# phaseType =1 unblinded only; =2 blinded only; =3 unblinded & blinded
#
setGeneric(name = ".verifyPhase",
           def = function(phase, ...) { 
               standardGeneric(".verifyPhase") 
             })

# anything not explicitly allowed is forbidden
setMethod(f = ".verifyPhase",
          signature = c(phase = "ANY"),
          definition = function(phase, ...) { 
              stop("inappropriate object provided for phase", call. = FALSE)
            })

setMethod(f = ".verifyPhase",
          signature = c(phase = "character"),
          definition = function(phase, ..., cutoff, data, subset, dObj) { 

              if (!is.numeric(x = cutoff)) {
                stop("cutoff must numeric [0,1]", call. = FALSE)
              }
              if (cutoff < 0.0 || cutoff > 1.0) {
                stop("cutoff must numeric [0,1]", call. = FALSE)
              }

              phase <- tolower(x = phase)

              # times of infection for blinded stage
              # limited to only those with the variant(s) under analysis

              grp <- subset & {data[,dObj$Gam] == 0L}
              infTimesBlinded <- sort(x = unique(x = data[grp, dObj$U]))
              nInfTimesB <- length(x = infTimesBlinded)

              # times of infection for unblinded stage
              # limited to only those with the variant(s) under analysis

              grp <- subset & {data[,dObj$Gam] == 1L}
              infTimesUnblinded <- sort(x = unique(x = data[grp, dObj$U]))
              nInfTimesU <- length(x = infTimesUnblinded)

              if (phase %in% c("ub", "bu")) {

                if (nInfTimesB == 0L && nInfTimesU == 0L) {
                  stop("no infections identified", call. = FALSE)
                }

                # both phases to be included in analysis -- make sure 
                # sufficient #s are present
                if (nInfTimesB < max(round(cutoff*{nInfTimesB + nInfTimesU}),1.0)) {

                  message("< ", round(x = cutoff*100, digits = 2),
                          "% (", nInfTimesB, "/", {nInfTimesB + nInfTimesU},
                          ") of infections occurred during the unblinded phase;\n",
                          "only the blinded phase is included in analysis")
                  infTimesBlinded <- NULL
                  nInfTimesB <- 0L

                  phaseType <- 1L

                } else if (nInfTimesU < max(round(cutoff*{nInfTimesB + nInfTimesU}),1.0)) {

                  message("< ", round(x = cutoff*100, digits = 2),
                          "% (", nInfTimesU, "/", {nInfTimesB + nInfTimesU},
                          ") of infections occurred during the unblinded phase;\n",
                          "only the blinded phase is included in analysis")
                  infTimesUnblinded <- NULL
                  nInfTimesU <- 0L

                  phaseType <- 2L

                } else {

                  message("both blinded and unblinded phases included in analysis")
                  phaseType <- 3L

                }
              } else if (phase == "u") {

                  if (nInfTimesU == 0L) {
                    stop("no infections occurred during unblinded phase", 
                         call. = FALSE)
                  }

                  infTimesBlinded <- NULL
                  nInfTimesB <- 0L
                  message("only the unblinded phase is included in analysis")
                  phaseType <- 1L

              } else if (phase == "b") {

                  if (nInfTimesB == 0L) {
                    stop("no infections occurred during blinded phase", 
                         call. = FALSE)
                  }

                  infTimesUnblinded <- NULL
                  nInfTimesU <- 0L
                  message("only the blinded phase is included in analysis")
                  phaseType <- 2L

              }

              return( list("phaseType" = phaseType,
                           "infTimesBlinded" = infTimesBlinded,
                           "infTimesUnblinded" = infTimesUnblinded) )
            })
