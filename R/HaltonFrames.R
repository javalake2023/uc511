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
  # validate the parm being validated i.e. is it a supported uc511 function parameter.
  if(!parm %in% base::c("n", "J", "bases", "shapefile", "panels", "panel_overlap", "randomStart",
                  "shp", "bb", "stratum", "nExtra", "quiet", "inclSeed",
                  "seeds",
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


#' @name HaltonFrame
#'
#' @title Generate a Halton Frame.
#'
#' @description A description of this useful function.
#'
#' @details This function was written by Phil Davies.
#'
#' @param n The number of points in the frame to generate.
#' @param J The number of grid cells. A list of 2 values. The default value is c(3, 2), we could also use c(5, 3).
#' @param bases Co-prime base for the Halton Sequence. The default value is c(2, 3).
#' @param shapefile A sf object. If the shapefile parameter is NULL then function
#' uc511::HaltonFrameBase is called directly.
#' @param panels A list of integers that define the size of each panel in a
#' non-overlapping panels design. The length of the list determines the number of
#' panels required. The sum of the integers in the panels parameter will determine
#' the total number of samples selected, n. The default value for panels is NULL,
#' this indicates that a non-overlapping panel design is not wanted.
#' @param panel_overlap A list of integers that define the overlap into the previous
#' panel. Is only used when the panels parameter is not NULL. The default value for
#' panel_overlap is NULL. The length of panel_overlap must be equal to the length
#' of panels. The first value is always forced to zero as the first panel never
#' overlaps any region.
#' @param randomStart Whether a spatially balanced sample will be randomly drawn from
#' the frame or not. Default value is FALSE.
#' @param seeds A list of 2 seeds, u1 and u2. If not specified, default is NULL.
#' @param stratum Name of column in shapefile that makes up the strata.
#'
#' @return A list containing the following variables:
#'         - halton_seq
#'         - halton_seq_div
#'         - Z
#'         - halton_frame
#'         - J
#'         - sample
#'         - pts.shp
#'         - bb
#'
#' @examples
#' \dontrun{
#' hf_ <- uc511::HaltonFrame()
#' }
#'
#' @export
HaltonFrame <- function(n = (bases[1]^J[1]) * (bases[2]^J[2]),
                        J = base::c(3, 2),
                        bases = base::c(2, 3),
                        shapefile = NULL,
                        panels = NULL,
                        panel_overlap = NULL,
                        randomStart = FALSE,
                        seeds = NULL,
                        stratum = NULL){

  # validate our parameters.
  uc511::validate_parameters("J", J)
  uc511::validate_parameters("bases", bases)
  uc511::validate_parameters("n", c(n))
  if (!is.null(seeds)){
    uc511::validate_parameters("seeds", seeds)
  } else {
    seeds <- base::floor(stats::runif(2, 0, 62208))
  }

  # validate panel design if we are using one.
  res <- ValidatePanelDesign(panels, panel_overlap, n)
  panel_design  <- res$panel_design
  number_panels <- res$number_panels
  panel_overlap <- res$panel_overlap
  n             <- res$n

  # panel_design and randomStart are mutually exclusive. stop if both TRUE.
  if(panel_design & randomStart){
    stop("uc511(HaltonFrame) Panel design and randomStart are mutually exclusive.")
  }

  # state how many samples user is looking for.
  msg <- "uc511(HaltonFrame) %s samples have been requested."
  msgs <- sprintf(msg, n)
  message(msgs)

  # validate the shapefile (if specified) has an associated CRS.
  crs <- NULL
  # Check if the shapefile has an associated CRS
  if (!is.null(shapefile)){
    if (is.null(sf::st_crs(shapefile))) {
      stop("uc511(HaltonFrame) Shapefile does not have an associated CRS.")
    } else {
      msg <- "uc511(HaltonFrame) Shapefile has an associated CRS."
      msgs <- sprintf(msg)
      message(msgs)
      crs <- sf::st_crs(shapefile)
    }
  } else {
    # shapefile is null, ie. not specified, so just call HaltonFrameBase
    hf_ <- uc511::HaltonFrameBase(J = c(J[1], J[2]), bases = bases, seeds = seeds)
    return(hf_)
  }

  # number of points currently in the area of interest.
  pts_in_intersection <- 0
  # initialise
  i <- 0

  while (pts_in_intersection <= n){
    # keep going until we have found the required number of points

    # create halton frame
    hf_ <- uc511::HaltonFrameBase(J = c(J[1]+i, J[2]+i), bases = bases, seeds = seeds)

    # process points returned.
    pts <- hf_$halton_frame
    pts <- cbind(seq(1, dim(pts)[1]), pts)

    # save returned seeds in case they have changed (would only change if initially NULL).
    seeds <- hf_$seeds

    #
    bb <- sf::st_as_sfc(sf::st_bbox(shapefile))
    cntrd <- sf::st_centroid(bb)
    bb.rot <- (bb - cntrd) * rot(0) + cntrd
    bb.new <- sf::st_as_sfc(sf::st_bbox(bb.rot))

    #
    base::attr(bb.new, "rotation") <- 0
    base::attr(bb.new, "centroid") <- sf::st_coordinates(cntrd)
    pts.shp <- uc511::rotate.scale.coords(coords = pts, bb = bb.new)

    # replace the CRS (or set), st_intersection needs both objects with the same CRS.
    sf::st_crs(shapefile) <- crs
    sf::st_crs(pts.shp) <- crs
    # make the assumption (that the attribute is constant throughout the geometry) explicit# make the assumption (that the attribute is constant throughout the geometry)
    #sf::st_agr(shp.ashb) <- "constant"
    #sf::st_agr(pts.shp) <- "constant"

    #tmp <- sf::st_as_sf(pts.shp)
    tmp <- sf::st_cast(pts.shp, "POINT")
    tmp <- sf::st_as_sf(tmp)
    tmp$ID <- seq(1, dim(pts)[1])
    diff_ <- sf::st_intersection(tmp, shapefile)

    # find number of points within our shapefile.
    pts_in_intersection <- length(sf::st_cast(sf::st_union(diff_), "POINT"))

    msg <- "uc511(HaltonFrame) points in intersection: %s."
    msgs <- sprintf(msg, pts_in_intersection)
    message(msgs)

    # expand the Halton frame.
    i <- i + 1
  }

  # display some statistics and return results.
  msg <- "uc511(HaltonFrame) %s samples found in %s iterations, using J1=%s and J2=%s."
  msgs <- sprintf(msg, pts_in_intersection, i, J[1]+i-1, J[2]+i-1)
  message(msgs)

  # are we performing a randomStart?
  if (randomStart){
    # need to select n samples from diff_ (rename this to sample).
    # and then replicate; generate random number and take new sample.
    message("uc511(HaltonFrame) randomStart.")
    diff_pts <- sf::st_cast(diff_, "POINT")
    df_sorted <- diff_pts[order(diff_pts$ID),]
    df_sorted <- sf::st_as_sf(df_sorted)
    df_sorted$uc511SeqID <- seq(1, length(df_sorted$ID))
    duplicated_pts <- base::rbind(df_sorted, df_sorted)
    random_start_point <- base::sample(1:length(duplicated_pts$ID), 1)
    sample_indices <- base::seq(random_start_point, (random_start_point + n) - 1, 1)
    random_start_sample <- duplicated_pts[sample_indices,]
    diff_ <- random_start_sample
    # randomStart mutually exclusive with panel_design.
    panel_design <- FALSE
  }

  # are we performing a panel_design? yes then go assign panelid's.
  if(panel_design){
    message("uc511(HaltonFrame) panel_design.")
    diff_pts <- sf::st_cast(diff_, "POINT")
    df_sorted <- diff_pts[order(diff_pts$ID), ]
    df_sorted$uc511SeqID <- seq(1, length(df_sorted$ID))
    diff_pts_sf <- sf::st_as_sf(df_sorted)
    res <- PanelDesignAssignPanelids(diff_pts_sf, panels, panel_overlap, panel_design, number_panels)
    diff_ <- res$sample
  }

  # if we are not performing a randomStart or a panel_design
  if (!randomStart & !panel_design){
    message("uc511(HaltonFrame) return n points.")
    # turn our sample into points.
    diff_pts <- sf::st_cast(diff_, "POINT")
    df_sorted <- diff_pts[order(diff_pts$ID), ]
    df_sorted$uc511SeqID <- seq(1, length(df_sorted$ID))
    # return everything from the intersection.
    diff_ <- df_sorted #[1:n,]
  }

  # Need to return cpprshs$pts, cpprshs$xklist, z and hf
  result <- base::list(halton_seq     = hf_$halton_seq,
                       halton_seq_div = hf_$halton_seq_div,
                       Z              = hf_$Z,
                       halton_frame   = hf_$halton_frame,
                       J              = c(J[1]+i-1, J[2]+i-1),
                       sample         = diff_,
                       hf.pts.shp     = pts.shp,
                       bb             = bb.new,
                       seeds          = seeds)
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
#' @param J The number of grid cells. A list of 2 values. The default value is c(3, 2), we could also use c(5, 3).
#' @param bases Co-prime base for the Halton Sequence. The default value is c(2, 3).
#' @param seeds The u1 and u2 seeds to use.
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
HaltonFrameBase <- function(n = (bases[1]^J[1]) * (bases[2]^J[2]),
                            J = base::c(3, 2),
                            bases = base::c(2, 3),
                            seeds = NULL){
  # validate our parameters.
  uc511::validate_parameters("J", J)
  uc511::validate_parameters("bases", bases)
  uc511::validate_parameters("n", base::c(n))
  if (!is.null(seeds)){
    uc511::validate_parameters("seeds", seeds)
  }

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
  if(is.null(seeds)){
    cpprshs <- uc511::cppBASpts(n = B2, bases = bases)
  } else {
    cpprshs <- uc511::cppBASpts(n = B2, bases = bases, seeds = seeds)
  }
  #
  z <- ((cpprshs$pts[1:B,] + cpprshs$pts[(B+1):B2,])/2)
  # x-dimension
  x_dim <- (base::floor(bases[1]^j1 * z[,1])/bases[1]^j1) + 0.5*(1/(bases[1]^j1))
  # y-dimension
  y_dim <- (base::floor(bases[2]^j2 * z[,2])/bases[2]^j2) + 0.5*(1/bases[2]^j2)
  # Halton Frame
  hf <- base::cbind(x_dim, y_dim)
  # Create an index number for each data point.
  halton_indx <- base::seq(1, B2/2)

  # Need to return cpprshs$pts, cpprshs$xklist, z, hf and seeds.
  result <- base::list(halton_seq     = cpprshs$pts,
                       halton_seq_div = cpprshs$xklist,
                       Z              = z,
                       halton_frame   = hf,
                       halton_indx    = halton_indx,
                       seeds          = cpprshs$seeds)
  return(result)
}
