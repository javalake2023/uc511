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
#' @param shapefile something
#' @param crs something
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
HaltonFrame <- function(n = (bases[1]^J[1]) * (bases[2]^J[2]),
                        J = c(3, 2),
                        bases = c(2, 3),
                        shapefile = NULL,
                        crs = NULL){
  # validate our parameters.
  uc511::validate_parameters("J", J)
  uc511::validate_parameters("bases", bases)
  uc511::validate_parameters("n", c(n))

  # number of points currently in the area of interest.
  pts_in_intersection <- 0
  # initialise
  i <- 0

  while (pts_in_intersection <= n){
    # keep going until we have found the required number of points

    # create halton frame
    hf_ <- uc511::HaltonFrameBase(J = c(J[1]+i, J[2]+i), bases = bases)

    # process points returned.
    pts <- hf_$halton_frame
    pts <- cbind(seq(1, dim(pts)[1]), pts)

    #
    bb <- sf::st_as_sfc(sf::st_bbox(shapefile))
    cntrd <- sf::st_centroid(bb)
    bb.rot <- (bb - cntrd) * rot(0) + cntrd
    bb.new <- sf::st_as_sfc(sf::st_bbox(bb.rot))

    #
    attr(bb.new, "rotation") <- 0
    attr(bb.new, "centroid") <- sf::st_coordinates(cntrd)
    pts.shp <- uc511::rotate.scale.coords(coords = pts, bb = bb.new)

    # replace the CRS (or set), st_intersection needs both objects with the same CRS.
    sf::st_crs(shapefile) <- crs
    sf::st_crs(pts.shp) <- crs
    # make the assumption (that the attribute is constant throughout the geometry) explicit# make the assumption (that the attribute is constant throughout the geometry)
    #sf::st_agr(shp.ashb) <- "constant"
    #sf::st_agr(pts.shp) <- "constant"
    diff_ <- sf::st_intersection(shapefile, pts.shp)
    # find number of points within our shapefile.
    pts_in_intersection <- lengths(diff_$geometry)
    i <- i + 1
  }

  # display some statistics and return results.
  msg <- "uc511(HaltonFrame) %s samples found in %s iterations, using J1=%s and J2=%s."
  msgs <- sprintf(msg, pts_in_intersection, i, J[1]+i, J[2]+i)
  message(msgs)

  # Need to return cpprshs$pts, cpprshs$xklist, z and hf
  result <- base::list(halton_seq = hf_$halton_seq,
                       halton_seq_div = hf_$halton_seq_div,
                       Z = hf_$Z,
                       halton_frame = hf_$halton_frame,
                       J = c(J[1]+i, J[2]+i),
                       diff_ = diff_,
                       pts.shp = pts.shp,
                       bb = bb.new)
  return(result)
}


#' @name HaltonFrameBase
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
#' hf_ <- uc511::HaltonFrameBase()
#' }
#'
#' @export
HaltonFrameBase <- function(n = (bases[1]^J[1]) * (bases[2]^J[2]), J = c(3, 2), bases = c(2, 3)){
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
    B <- (base::floor(n / B) + 1) * B
  }
  # double the number of points
  B2 <- 2 * B
  # compute B2 halton points
  #cpprshs <- uc511::cppRSHalton_br(n = B2, bases = bases)
  cpprshs <- uc511::cppBASpts(n = B2, bases = bases)
  #
  z <- ((cpprshs$pts[1:B,] + cpprshs$pts[(B+1):B2,])/2)
  # x-dimension
  x_dim <- (base::floor(bases[1]^j1 * z[,1])/bases[1]^j1) + 0.5*(1/(bases[1]^j1))
  # y-dimension
  y_dim <- (base::floor(bases[2]^j2 * z[,2])/bases[2]^j2) + 0.5*(1/bases[2]^j2)
  # Halton Frame
  hf <- base::cbind(x_dim, y_dim)

  # Need to return cpprshs$pts, cpprshs$xklist, z and hf
  result <- base::list(halton_seq = cpprshs$pts,
                       halton_seq_div = cpprshs$xklist,
                       Z = z,
                       halton_frame = hf)
  return(result)
}
