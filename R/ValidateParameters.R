# ValidateParameters.R

#' @name validate_parameters
#'
#' @title Validate uc511 function parameters.
#'
#' @description This function is used to validate parameters passed to
#' all available uc511 functions.
#'
#' @details This function was written by Phil Davies.
#'
#' @param parm The parameter to be validated.
#' @param parm_value The value of the parameter to be validated. Must be defined as a list.
#'
#' @return Always returns TRUE indicating that the parameter was parsed successfully. If
#' a parameter fails validation further execution is terminated using the STOP function.
#'
#' @export
validate_parameters <- function(parm, parm_value){
  # validate the parm being validated i.e. is it a supported uc511 function parameter.
  if(!parm %in% base::c("n", "J", "bases", "shapefile", "panels", "panel_overlap",
                        "randomStart", "N",
                        "shp", "bb", "stratum", "nExtra", "quiet", "inclSeed",
                        "seeds", "boundingbox",
                        "seed", "total_rows", "sample_size",
                        "testparm",
                        "panelid",
                        "hipPopulation", "hipN", "hipIterations")){
    stop(base::c("uc511(validate_parameters) specified parameter ", parm, " is not currently supported."))
  }

  # make sure parameter is a list (or vector)
  if(!is.vector(parm_value)){
    stop(base::c("uc511(validate_parameters) Parameter ", parm, " must be a list or vector."))
  }
  # check if all values in the list are numeric
  if(!all(sapply(parm_value, is.numeric))){
    stop(base::c("uc511(validate_parameters) Parameter ", parm, " must contain all numeric values."))
  }
  # check if list is of length 2
  if(parm == "J" & length(parm_value) != 2){
    stop("uc511(validate_parameters) Parameter J must be a list of length 2.")
  }
  # check if values for SRS parameters are greater than zero. check if panelid is greaer than zero.
  if(parm %in% base::c("seed", "total_rows", "sample_size", "panelid")) {
    if (parm_value <= 0){
      stop(base::c("uc511(validate_parameters) Parameter ", parm, " must have a value greater than zero."))
    }
  }
  # validate HIP parameters
  if(parm == "hipIterations"){
    if(parm_value < 2){
      stop(base::c("uc511(validate_parameters) Parameter ", parm, " values less than two are not supported."))
    }
    if(parm_value > 13){
      stop(base::c("uc511(validate_parameters) Parameter ", parm, " values greater than 13 are not supported."))
    }
  }
  return(TRUE)
}
