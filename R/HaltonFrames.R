# HaltonFrames.R

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
  # make sure parameter is a list (or vector)
  if(!is.vector(parm_value)){
    stop(c("uc511(validate_parameters) Parameter ", parm, " must be a list or vector."))
  }
  # check if all values in the list are numeric
  if(!all(sapply(parm_value, is.numeric))){
    stop(c("uc511(validate_parameters) Parameter ", parm, " must contain all numeric values."))
  }
  # check if list is of length 2
  if(parm == "J" & length(parm_value) != 2){
    stop("uc511(validate_parameters) Parameter J must be a list of length 2.")
  }
  # check if values for SRS parameters are greater than zero.
  if(parm %in% c("seed", "total_rows", "sample_size")) {
    if (parm_value <= 0){
      stop(c("uc511(validate_parameters) Parameter ", parm, " must have a value greater than zero."))
    }
  }
  # validate HIP parameters
  if(parm == "hipIterations"){
    if(parm_value < 2){
      stop(c("uc511(validate_parameters) Parameter ", parm, " values less than two are not supported."))
    }
    if(parm_value > 13){
      stop(c("uc511(validate_parameters) Parameter ", parm, " values greater than 13 are not supported."))
    }
  }
  return(TRUE)
}


#' @name HaltonFrame
#'
#' @title Generate a Halton Frame.
#'
#' @description A description of this useful function.
#'
#' @details This function was written by Phil Davies.
#'
#' @param n The number of points in the frame to generate.
#' @param J A list of 2 values. The default value is c(3, 2), we could also use c(5, 3).
#' @param bases Co-prime base for the Halton Sequence. The default value is c(2, 3).
#'
#' @return A list containing the following four variables:
#' halton_seq -
#' halton_seq_div -
#' Z -
#' halton_frame -
#'
#' @examples
#' \dontrun{
#' hf_ <- uc511::HaltonFrame()
#' }
#'
#' @export
HaltonFrame <- function(n = (bases[1]^J[1]) * (bases[2]^J[2]), J = c(3, 2), bases = c(2, 3)){
  # validate our parameters.
  uc511::validate_parameters("J", J)
  uc511::validate_parameters("bases", bases)
  uc511::validate_parameters("n", c(n))

  #
  j1 <- J[1]; j2 <- J[2]
  # calculate B
  B <- (bases[1]^j1) * (bases[2]^j2)
  # check how many points the caller wants.
  if(n > B) {
    B <- (floor(n / B) + 1) * B
  }
  # double the number of points
  B2 <- 2 * B
  # compute B2 halton points
  #cpprshs <- uc511::cppRSHalton_br(n = B2, bases = bases)
  cpprshs <- uc511::cppBASpts(n = B2, bases = bases)
  #
  z <- ((cpprshs$pts[1:B,] + cpprshs$pts[(B+1):B2,])/2)
  # x-dimension
  x_dim <- (floor(bases[1]^j1 * z[,1])/bases[1]^j1) + 0.5*(1/(bases[1]^j1))
  # y-dimension
  y_dim <- (floor(bases[2]^j2 * z[,2])/bases[2]^j2) + 0.5*(1/bases[2]^j2)
  # Halton Frame
  hf <- cbind(x_dim, y_dim)

  # Need to return cpprshs$pts, cpprshs$xklist, z and hf
  result <- list(halton_seq = cpprshs$pts,
                 halton_seq_div = cpprshs$xklist,
                 Z = z,
                 halton_frame = hf)
  return(result)
}
