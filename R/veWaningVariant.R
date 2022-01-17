#' Estimation of Vaccine Efficacy Over Time - Variant Aware
#'
#' Implements methods for inference on potential waning of vaccine
#'    efficacy and for estimation of vaccine efficacy at a user-specified time
#'    after vaccination based on data from a randomized, double-blind,
#'    placebo-controlled vaccine trial in which participants may be unblinded
#'    and placebo subjects may be crossed over to the study vaccine.  The
#'    method allows for variant specification and adjustment for possible 
#'    confounding via inverse
#'    probability weighting through specification of models for the trial
#'    entry process, unblinding mechanisms, censoring process, and the 
#'    probability an unblinded placebo participant accepts study vaccine. 
#'
#' Analysis data must include the following information:
#'   \describe{
#'   \item{E}{The study entry time.}
#'   \item{A}{Binary treatment assignment, where A=0 indicates placebo; 
#'            A = 1 otherwise.}
#'   \item{U}{The minimum of: time to infection, time to refusing study 
#'      vaccine after unblinding (placebo participants only), or time to censoring 
#'      (loss to followup or administrative censoring due to analysis at time L).}
#'   \item{Delta}{Infection-variant category indicator, where: 
#'     -1 indicates a placebo participant was unblinded and refused study vaccine;
#'     0 indicates censoring; and
#'     = ve indicates infection of variant ve (ve = 1, ..., nV).}
#'   \item{R}{The minimum of: time to unblinding, time to infection, or
#'     time to censoring.}
#'   \item{Gamma}{Indicator of R dynamic, where 1 indicates that R is the
#'     time to unblinding and 0 indicates that R is either the time to infection
#'     or the time to censoring.}
#'   \item{Psi}{Indicator of whether a subject received study vaccine, where 1
#'     indicates that a participant was assigned to vaccine or was assigned to 
#'     placebo is unblinded and decides to get the study vaccine; 0 otherwise.}
#' }
#' Data can include baseline covariates. Methods for time-dependent covariates
#'  are not currently available.
#'
#' If input lag is provided as a numeric vector or as a column of the data,
#'   lag values should be set to Inf or NA for participants that do not reach 
#'   full efficacy.
#'
#' There are 3 possible weighting options, the selection of which is determined
#'  by the combination of model inputs. 
#'   \describe{
#'     \item{No weighting: }{No models are provided as input, i.e., 
#'           inputs modelEntry, modelUnblind0, modelUnblind1, modelPsi, 
#'           modelCensor0, and modelCensor1 are not provided or are NULL.}
#'     \item{Weighting depends only on the censoring dynamic: }{Models
#'           modelCensor0 and modelCensor1 must be provided and models 
#'           modelEntry, modelUnblind0,  modelUnblind1, modelPsi must be NULL.}
#'     \item{Weighting depends on the entry, unblinding, and censoring 
#'           dynamics: }{All models must be provided.}
#'  }
#'
#' The returned S3 object has 6 attributes needed for post-processing
#'  tools ve() and plot(). Specifically, "gFunc" is an integer object
#'  specifying the model selected for the infection rate (based on input gFunc);
#'  "v", the knots or cut-offs to be used by gFunc (input v);
#'  "maxTime", the maximum time since full efficacy included in the analysis; 
#'  "phaseType", =1 only unblinded phase, =2, only blinded phase, =3 both phases,
#'  and "wgtType", =0 no weighting, =1 censor weighting, =2 full weighting.
#'
#' @param data A data.frame object containing all relevant data. 
#'   Records with missing data will be removed prior to initiating the analysis. 
#'
#' @param L A numeric object. The analysis time.
#'
#' @param ... Ignored. Used only to require named inputs.
#'
#' @param phase A character object. The phase(s) to include in the analysis.
#'   Default is ="ub", indicating that both blinded and unblinded phases will be
#'   used to estimate theta values if possible. If ="b", only the blinded
#'   phase will be used. If ="u", only the unblinded phase will be used.
#'
#' @param cutoff A numeric object. The minimum proportion of infections that
#'   must occur during a phase to be included in the analysis. The default
#'   is 0.1 (10 %). Specifically -- if < 10% of the infections occur during the 
#'   blinded (unblinded) phase, only the unblinded (blinded) phase will be 
#'   included in the analysis.
#'
#' @param modelEntry A formula object or NULL. The coxph model for entry times.
#'   The LHS is set as the appropriate Surv() object internally. If a LHS
#'   is provided, it is ignored. If NULL, inputs modelPsi, modelUnblind0,
#'   and modelUnblind1 must also be NULL. See Details for further information.
#'
#' @param modelUnblind0 A formula object or NULL. The coxph model for 
#'   unblinding/crossover for placebo (A=0) participants. If NULL, inputs 
#'   modelEntry, modelPsi, and modelUnblind1 must also be NULL. See Details 
#'   for further information.
#'
#' @param modelUnblind1 A formula object or NULL. The coxph model for unblinding
#'   for vaccinated (A=1) participants. If NULL, inputs modelEntry, modelPsi, 
#'   and modelUnblind0 must also be NULL. See Details for further information.
#'
#' @param modelPsi A formula object or NULL. The logistic model for the
#'   probability that a placebo participant (A = 0) is unblinded at R
#'   (Gamma = 1) and agrees to take the study vaccine (Psi = 1).
#'   If a LHS is provided, it is ignored. If NULL, inputs modelEntry, 
#'   modelUnblind0, and modelUnblind1 must also be NULL. See Details for 
#'   further information.
#'
#' @param modelCensor0 A formula object or NULL. The coxph model for censoring
#'   for placebo (A=0) participants. The LHS is set as the appropriate Surv() 
#'   object internally. If a LHS is provided, it is ignored. If NULL, all other 
#'   models must be NULL. See Details for further information.
#'
#' @param modelCensor1 A formula object or NULL. The coxph model for censoring
#'   for vaccinated (A=1) participants. The LHS is set as the appropriate Surv() 
#'   object internally. If a LHS is provided, it is ignored. If NULL, all other 
#'   models must be NULL. See Details for further information.
#'
#' @param gFunc A character vector object. The model of infection rates.
#'   Must be one of {'lin', 'piece', 'splin', 'spcub'} for the linear,
#'   piecewise constant, linear spline, and cubic spline models, respectively.
#'
#' @param variant An integer object. The variant for the analysis. If 0,
#'   all variants are included in the analysis.
#'
#' @param lag A scalar numeric, numeric vector, or character object. The lag 
#'   time(s) between the initial vaccine dose and full efficacy. If a scalar,
#'   the provided lag time applies to all participants. If a numeric vector,
#'   the vector contains the individual specific lag time for each participant
#'   (see details for further information).
#'   If character, the column header of the data containing the lag times.
#'   The default value is a scalar value of 6 weeks (42 days) -- NOTE this 
#'   assumes that the data are on the scale of weeks. 
#'
#' @param v A numeric vector object. 
#'   The knots or cut-offs to be used for input gFunc.
#'   If gFunc = 'lin', this input is ignored. For 'splin' and 'spcub', the
#'   knots of the spline on (0,L). For 'piece', the cut-offs on 
#'   (0,L). Note that this input should not include the extremes 0 and L.
#'
#' @param minWgt A numeric object. If not NULL, the minimum non-zero value a 
#'   weight can have, i.e., weight = max(minWgt, weight). If NULL, no
#'   lower truncation of weights is performed.
#'
#' @param maxWgt A numeric object. If not NULL, the maximum value a 
#'   weight can have, i.e., weight = min(maxWgt, weight). If NULL, no
#'   upper truncation of weights is performed.
#'
#' @param txName A character object. The header of the column of data 
#'   containing the treatment variable. Default value is 'A'.
#'   Treatment must be coded as 0/1, where 1 indicates that participant
#'   was vaccinated; 0 otherwise.
#'
#' @param U A character object. The header of the column of data 
#'   containing the minimum of time to infection, time to refusing study 
#'   vaccine after unblinding (placebo participants only), or time to censoring
#'   (due to loss to follow up or administrative censoring). 
#'
#' @param entryTime A character object. The header of the column of data 
#'   containing the time of entry into the study on the scale of the
#'   calendar time. Default value is 'E'. 
#'
#' @param R A character object. The header of the column of data 
#'   containing the minimum of: time to unblinding, time to infection, or
#'   time to censoring. 
#'
#' @param Gamma A character object. The header of the column of data 
#'   containing the category for the R dynamic. Default value is 'Gam'.
#'   Data must be 0/1, where 1 indicates that R is the time to unblinding;
#'   0 indicates that R is the infection time or the censoring time. 
#'
#' @param Psi A character object. The header of the column of data 
#'   containing the indicator of whether a participant received study 
#'   vaccine, where 1 indicates that a participant was assigned to placebo is 
#'   unblinded and decides to get the study vaccine or that a participant
#'   was assigned to vaccine; 0 otherwise.  Default value is 'Psi'.
#'
#' @param Delta A character object. The header of the column of data 
#'   containing the infection-variant category indicator. 
#'
#' @returns A an S3 object of class "VEwaningVariant", which comprises a list
#'   object containing 
#'   \item{theta}{A vector object containing the estimated theta parameters.}
#'   \item{cov}{The covariance estimated using the sandwich estimator.}
#'   \item{SE}{The standard error estimated using the sandwich estimator.}
#'   and attributes  "gFunc", "maxTime", "v", "phaseType", and "wgtType", 
#'   which store
#'   details of the original analysis that are required for post-processing
#'   convenience functions ve() and plot(). See details for further
#'   information.
#'
#' @export
#' @examples
#' data(variantData)
#'
#' set.seed(1234)
#'
#' ind <- sample(1:nrow(variantData), 2000)
#' # NOTE: This sample size is chosen for example only -- larger data sets
#' # should be used.
#'
#' # no weighting -- variant 1 infection only
#' 
#' res_noWgt <- veWaningVariant(data = variantData[ind,], 
#'                              L = 52.0,  
#'                              variant = 1L,
#'                              gFunc = 'piece', 
#'                              v = c(5.0,10.0))
#'
#' # censoring only weighting -- variant 1 infection only
#' 
#' res_cens <- veWaningVariant(data = variantData[ind,], 
#'                             L = 52.0,  
#'                             variant = 1L,
#'                             modelCensor0 = ~ X1+X2, 
#'                             modelCensor1 = ~ X1+X2, 
#'                             gFunc = 'piece', 
#'                             v = c(5.0,10.0))
#'
#' # full weighting -- variant 1 infection only
#' 
#' \dontrun{res_full <- veWaningVariant(data = variantData[ind,], 
#'                             L = 52.0,  
#'                             variant = 1L,
#'                             modelEntry = ~ X1,
#'                             modelUnblind0 = ~X1+X2,
#'                             modelUnblind1 = ~X2,
#'                             modelPsi = ~X1*X2,
#'                             modelCensor0 = ~ X1+X2, 
#'                             modelCensor1 = ~ X1+X2, 
#'                             gFunc = 'piece', 
#'                             v = c(5.0,10.0))}
#'
#' @import survival
#' @import Rcpp
#' @import methods
#' @include verifyInputs.R verifyModels.R verifyPhase.R 
#' @include fitERPsi.R fitCensor.R esttheta.R
#'
veWaningVariant <- function(data, 
                            L, 
                            ..., 
                            phase = "ub",
                            cutoff = 0.1,
                            lag = 6.0,
                            modelEntry = NULL,
                            modelUnblind0 = NULL,
                            modelUnblind1 = NULL,
                            modelPsi = NULL,
                            modelCensor0 = NULL,
                            modelCensor1 = NULL,
                            gFunc = NULL,
                            variant = 0L,
                            v = NULL,
                            minWgt = NULL,
                            maxWgt = NULL,
                            txName = "A",
                            U = "U", 
                            entryTime = "E", 
                            Gamma = "Gam", 
                            R = "R",
                            Psi = "Psi",
                            Delta = "Delta") {

  # data must be a data.frame or matrix with column headers
  if (is.matrix(x = data)) {
    if (is.null(x = colnames(x = data))) {
      stop("data must be a data.frame object", call. = FALSE)
    }
    data <- as.data.frame(x = data)
  }
  if (!is.data.frame(x = data)) {
      stop("data must be a data.frame object", call. = FALSE)
  }

  # ensure the provided column names are characters
  if (any(!is.character(x = txName) ||
          !is.character(x = U) ||
          !is.character(x = entryTime) ||
          !is.character(x = Gamma) ||
          !is.character(x = R) ||
          !is.character(x = Psi) ||
          !is.character(x = Delta))) {
    stop("{txName, U, entryTime, Gamma, R, Psi, Delta} must be characters",
         call. = FALSE)
  }

  # list to store the column names of relevant data
  dObj <- list("Gam" = Gamma, 
               "R" = R, 
               "E" = entryTime,
               "U" = U,
               "A" = txName,
               "Psi" = Psi,
               "Delta" = Delta)

  # generate random name for lag variable if not provided
  if (is.character(x = lag)) {
    dObj$lag <- lag
  } else {
    while (TRUE) {
      lagName = paste(sample(x = LETTERS, size = 10, replace = TRUE), 
                      collapse = "")
      if (lagName %in% colnames(x = data)) next
      break
    }
    dObj$lag <- lagName
  }

  if (is.numeric(x = lag)) {
    # This catches both integer and numeric objects
    if ({length(x = lag) == 1L} || {length(x = lag) == nrow(x = data)}) {
      data[,dObj$lag] <- lag
    } else {
      stop("when provided as a numeric, lag must be a scalar ",
           "or a vector of length nrow(data)", 
           call. = FALSE)
    }
  } else if (!is.character(x = lag)) {
    stop("lag must be numeric scalar, numeric vector, or character object", 
         call. = FALSE)
  }

  # extracting baseline covariates included in models
  # and ensuring appropriate combinations are provided
  # for the weighting
  modelInfo <- .verifyModels(entry = modelEntry,
                             unblind0 = modelUnblind0,
                             unblind1 = modelUnblind1,
                             psi = modelPsi,
                             censor0 = modelCensor0,
                             censor1 = modelCensor1)

  # verify data and g-function specifications
  inputs <- .verifyInputs(data = data, 
                          dObj = dObj,
                          covs = modelInfo$covs,
                          gFunc = gFunc,  
                          variant = variant,  
                          v = v,  
                          L = L)

  data <- inputs$data
  gFunc <- inputs$gFuncObj

  # verify phases to be included in analysis and extract infection times
  phaseInfo <- .verifyPhase(phase = phase, 
                            cutoff = cutoff,
                            data = data, 
                            subset = gFunc$group, 
                            dObj = dObj)

  fitERPsi <- .fitERPsi(data = data, 
                        entry = modelEntry,  
                        unblind0 = modelUnblind0,  
                        unblind1 = modelUnblind1,  
                        Psi = modelPsi,  
                        dObj = dObj,
                        infTimesBlinded = phaseInfo$infTimesBlinded,
                        wgtType = modelInfo$wgtType)

  fitCensor <- .fitCensor(data = data, 
                          censor0 = modelCensor0, 
                          censor1 = modelCensor1, 
                          dObj = dObj, 
                          wgtType = modelInfo$wgtType, 
                          infTimesBlinded = phaseInfo$infTimesBlinded, 
                          infTimesUnblinded = phaseInfo$infTimesUnblinded, 
                          L = L)

  if (is.null(x = minWgt)) minWgt <- 0.0
  if (is.null(x = maxWgt)) maxWgt <- 1e8
  if (minWgt < 0.0) stop("minWgt cannot be negative", call. = FALSE)
  if (maxWgt < 0.0) stop("maxWgt cannot be negative", call. = FALSE)
  if (maxWgt < minWgt) stop("maxWgt cannot be < minWgt", call. = FALSE)

  # get theta estimates, covariance matrix, and standard errors
  out <- esttheta(data = data,
                  timesB = phaseInfo$infTimesBlinded,
                  timesU = phaseInfo$infTimesUnblinded,
                  fitERPsi = fitERPsi,
                  censorFit = fitCensor,
                  gFunc = gFunc,
                  dObj = dObj,
                  minWgt = minWgt,
                  maxWgt = maxWgt,
                  wgtType = modelInfo$wgtType,
                  phaseType = phaseInfo$phaseType)

  # any post-processing tools must respect the upper limit of the training
  #   set
  if (phaseInfo$phaseType != 2L) {
    # for vaccine participants, tau is measured from entry time
    # for placebo participants that accepted vaccine after unblinding,
    #   tau is measured from the unblinding time
    allTau <- {data[,dObj$A] == 1L}*{L - data[,dObj$E]} +
              {data[,dObj$A] == 0L}*{data[,dObj$Psi] == 1L}*{L - data[,dObj$R]}
  } else {
    # for blinded phase only, maxTau is the maximum time to unblinding
    allTau <- data[,dObj$R] * {data[,dObj$Gam] == 1L}
  }

  # maxTime is used in post-processing as the upper-bound of the range of
  # allowed times. Switching to "time since full efficacy" or (tau - lag)
  # means this value needs to be shifted down.
  maxTime <- max(allTau - data[,dObj$lag])
  if (maxTime < 0.0) {
    message("no participant has reached full efficacy")
    maxTime <- NULL
  }

  attr(out, "gFunc") <- gFunc$gFunc
  attr(out, "maxTime") <- maxTime
  attr(out, 'v') <- gFunc$v
  attr(out, 'phaseType') <- phaseInfo$phaseType
  attr(out, 'wgtType') <- modelInfo$wgtType

  class(out) <- "VEwaningVariant"

  return( out )
    
}
