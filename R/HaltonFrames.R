# HaltonFrames.R

#' @name HaltonFrame
#'
#' @title Generate a Halton Frame.
#'
#' @description Halton frames discretize an areal resource into a spatially ordered grid,
#' where samples of consecutive frame points are spatially balanced. To generate Halton Frames,
#' uc511 requires a study region \code{shapefile} and the region’s \code{bounding box}.
#'
#' @details This function was written by Phil Davies.
#'
#' @param N The number of points in the frame to generate.
#' @param J The number of grid cells. A list of 2 values. The default value is c(3, 2).
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
#' @param seeds A vector of 2 seeds, u1 and u2. If not specified, the default is NULL.
#' @param stratum Name of column in shapefile that makes up the strata.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return A list containing five variables:
#'
#' \itemize{
#' \item \code{J} The number of grid cells. A list of 2 values that were used to generate this
#' Halton grid and frame.
#' \item \code{hg.pts.shp} Halton grid over the bounding box and study area.
#' \item \code{hf.pts.shp} Halton frame, the sample points within the study area.
#' \item \code{bb} The bounding box.
#' \item \code{seeds} The u1 and u2 seeds used to generate the sample.
#' }
#'
#' The sample points in \code{hf.pts.shp} are returned in the form of a simple feature
#' collection of POINT objects. As well as having the features from the original \code{shapefile},
#' the following new attributes have been added:
#'
#' \itemize{
#'   \item \code{uc511SeqID}: A unique identifier for every sample point.
#'   \item \code{ID}: A unique identifier, the Halton frame point order.
#' }
#'
#' @examples
#' \dontrun{
#' hf_ <- uc511::HaltonFrame()
#' }
#'
#' @export
HaltonFrame <- function(N = 1,
                        J = base::c(3, 2),
                        bases = base::c(2, 3),
                        shapefile = NULL,
                        panels = NULL,
                        panel_overlap = NULL,
                        seeds = NULL,
                        stratum = NULL,
                        verbose = FALSE){

  # initialize variables.
  # defaults, before we see what the user wants.
  wantHaltonGrid <- FALSE
  wantHaltonFrame <- FALSE
  hf_stratification <- FALSE

  if(is.null(N)){
    wantHaltonGrid <- TRUE
    wantHaltonFrame <- FALSE
    if(verbose){
      base::message("uc511(HaltonFrame) Request for a Halton Grid.")
    }
  }

  # validate our parameters.
  uc511::validate_parameters("J", J)
  uc511::validate_parameters("bases", bases)
  if (!is.null(N)){
    uc511::validate_parameters("N", base::c(N))
  }

  if (is.null(seeds)){
    seeds <- uc511::generateUVector()
  } else {
    uc511::validate_parameters("seeds", seeds)
  }

  # validate panel design if we are using one.
  res <- ValidatePanelDesign(panels, panel_overlap, N)
  panel_design  <- res$panel_design
  number_panels <- res$number_panels
  panel_overlap <- res$panel_overlap
  N             <- res$n

  # if both not NULL then we want stratification.
  if(!base::is.null(base::names(N)) & !base::is.null(stratum)){
    hf_stratification <- TRUE
    wantHaltonFrame <- TRUE
    strata.levels <- base::names(N)
    if(verbose){
      msg <- "uc511(HaltonFrame) Stratification request for the following strata: %s.\n"
      msgs <- sprintf(msg, strata.levels)
      base::message(msgs)
    }
  } else {
    if(N >= 1){
      wantHaltonGrid <- TRUE
      wantHaltonFrame <- TRUE
      if(verbose){
        # state how many samples user is looking for.
        msg <- "uc511(HaltonFrame) Request for %s samples from a Halton Frame."
        msgs <- sprintf(msg, N)
        base::message(msgs)
      }
    }
  }

  # ensure the shapefile (if specified) has an associated CRS.
  crs <- NULL
  if (!is.null(shapefile)){
    if (is.null(sf::st_crs(shapefile))) {
      stop("uc511(HaltonFrame) Shapefile does not have an associated CRS.")
    } else {
      if(verbose){
        msg <- "uc511(HaltonFrame) Shapefile has an associated CRS."
        msgs <- sprintf(msg)
        base::message(msgs)
      }
      crs <- sf::st_crs(shapefile)
    }
  } else {
    # shapefile is null, ie. not specified, so just call HaltonFrameBase
    hf_ <- uc511::HaltonFrameBase(J = base::c(J[1], J[2]), bases = bases, seeds = seeds)
    return(hf_)
  }

  # number of points currently in the area of interest.
  pts_in_intersection <- 0
  # initialise
  i <- 0

  if(wantHaltonGrid & !wantHaltonFrame){
    # go get halton grid.
    result <- getHaltonFrame(shapefile, J, i, bases, seeds, crs)
    hf_ <- result$hf_
    diff_ <- result$sample
    pts.shp <- result$pts.shp
    bb.new <- result$bb.new
    seeds <- result$seeds

  } else {

    # need to check hf_stratification before we run the while loop.
    # run the while loop in new function so we can loop for each strata level for the desired N.
    # save the returned points using rbind.
    if(hf_stratification){
      # perform stratification
      smp <- NULL

      for(h in 1:base::length(N)){
        if(verbose){
          msg <- "uc511(BAS) Stratum: %s."
          msgs <- base::sprintf(msg, strata.levels[h])
          base::message(msgs)
        }

        h.indx <- base::which(shapefile[, stratum, drop = TRUE] == strata.levels[h])
        shp.stratum <- shapefile[h.indx,]

        # find the first point in the study region (picked at random)
        first.pt <- uc511::findFirstStudyRegionPoint(shapefile = shp.stratum, seeds = seeds)
        #
        k <- first.pt$k
        # generate seeds for the remaining points we need.
        seedshift <- base::c(first.pt$seeds[1] + k - 1, first.pt$seeds[2] + k - 1)

        result <- uc511::getHaltonPointsFromExpandableGrid(shapefile = shp.stratum,
                                                           N = N[h],
                                                           J = J,
                                                           bases = bases,
                                                           seeds = seedshift,
                                                           crs = crs,
                                                           verbose = verbose)

        seedshift <- result$seed
        diff_pts <- result$diff_
        df_sorted <- diff_pts[base::order(diff_pts$ID), ]
        #
        ret_sample <- rbind(first.pt$first.pt, df_sorted)
        sorted_samp <- ret_sample
        sorted_samp$uc511SeqID <- base::seq(1, base::length(sorted_samp$ID))
        diff_ <- sorted_samp[1:N[h],]
        # return original seeds.
        seeds <- first.pt$seeds
        #
        #df_sorted$uc511SeqID <- base::seq(1, base::length(df_sorted$ID))
        # return N[h] from the intersection for current stratum.
        #diff_ <- df_sorted[1:N[h],]
        smp <- base::rbind(smp, diff_)

      } # end for h
      # load variables for return to caller.
      i <- result$i
      diff_ <- smp
      pts.shp <- result$pts.shp
      bb.new <- result$bb.new
      #seeds <- result$seeds
      seeds <- first.pt$seeds

    } else {
      # stratification not required - use entire study area.

      # find the first point in the study region (picked at random)
      first.pt <- uc511::findFirstStudyRegionPoint(shapefile = shapefile, seeds = seeds)
      #
      k <- first.pt$k
      # generate seeds for the remaining points we need.
      seedshift <- base::c(first.pt$seeds[1] + k - 1, first.pt$seeds[2] + k - 1)
      # go get the Halton points.
      result <- uc511::getHaltonPointsFromExpandableGrid(shapefile = shapefile,
                                                         N = N,
                                                         J = J,
                                                         bases = bases,
                                                         seeds = seedshift,
                                                         crs = crs,
                                                         verbose = verbose)
      i         <- result$i
      diff_     <- result$diff_
      pts.shp   <- result$pts.shp
      bb.new    <- result$bb.new
      seedshift <- result$seeds         # was seeds until first.pt code was added.
      #pts_in_intersection <- result$pts_in_intersection
    }

  } # end if (wantHaltonGrid & !wantHaltonFrame

  # are we performing a panel_design? yes then go assign panelid's.
  if(panel_design){
    if(verbose){
      base::message("uc511(HaltonFrame) panel_design.")
    }
    diff_pts <- sf::st_cast(diff_, "POINT")
    df_sorted <- diff_pts[base::order(diff_pts$ID), ]
    df_sorted$uc511SeqID <- base::seq(1, base::length(df_sorted$ID))
    diff_pts_sf <- sf::st_as_sf(df_sorted)
    res <- PanelDesignAssignPanelids(diff_pts_sf, panels, panel_overlap, panel_design, number_panels)
    diff_ <- res$sample
  }

  # if we are not performing a randomStart or a panel_design or stratification
  if (!panel_design & wantHaltonFrame & !hf_stratification){
    if(verbose){
      message("uc511(HaltonFrame) Return N sample points.")
    }
    # turn our sample into points.
    diff_pts <- sf::st_cast(diff_, "POINT")
    df_sorted <- diff_pts[base::order(diff_pts$ID), ]
    #df_sorted$uc511SeqID <- base::seq(1, base::length(df_sorted$ID))   # move this down to after the rbind.
    # return everything from the intersection.
    diff_ <- df_sorted #[1:n,]

    #browser()
    # handle first.pt
    #f.pt <- sf::st_cast(first.pt$first.pt, "POINT")
    #f.pt <- sf::st_as_sf(f.pt)
    #zzz <- sf::st_as_sf(base::data.frame(SiteID = f.pt$ID, f.pt$x))
    ret_sample <- rbind(first.pt$first.pt, diff_)
    #sorted_samp <- ret_sample[base::order(ret_sample$SiteID), ]
    sorted_samp <- ret_sample
    sorted_samp$uc511SeqID <- base::seq(1, base::length(sorted_samp$ID))
    diff_ <- sorted_samp[1:n,]
    # return original seeds.
    seeds <- first.pt$seeds
  }

  # Need to return cpprshs$pts, cpprshs$xklist, z and hf
  result <- base::list(J              = c(J[1]+i-1, J[2]+i-1),
                       hf.pts.shp     = diff_,   # Halton Frame
                       hg.pts.shp     = pts.shp, # Halton Grid
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
    message("HaltonFrameBase seeds=", seeds)
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


#' @name getSample
#'
#' @title Generate a Halton Frame.
#'
#' @description A description of this useful function.
#'
#' @details This function was written by Phil Davies.
#'
#' @param shapefile A MULTIPOINT or POINT object from where to take the sample.
#' @param n The number of sample points to return.
#' @param randomStart Whether a spatially balanced sample will be randomly drawn from
#' the frame or not. Default value is FALSE.
#'
#' @return A list containing the following variable:
#'
#' \itemize{
#' \item \code{sample} The sample from the shapefile POINTS.
#' }
#'
#' @examples
#' \dontrun{
#' hf_ <- uc511::getSample()
#' }
#'
#' @export
# Drop the random sample functionality from HaltonFrame (sorry!) and have a 'getSample'
# function (like getPanel). Let 'Frame' be the output from HaltonFrame, then
#  - getSample(Frame, n): output the first n points from Frame
#  - getSample(Frame, n, randomStart = TRUE) output n points (in the order of the frame)
# from a random starting point in Frame.
getSample <- function(shapefile, n, randomStart = FALSE){
  # can accept input of either a MULTIPOINT object or a POINT object.
  # either $sample or $hf.pts.shp from the HaltonFrame function.
  uc511::validate_parameters("n", base::c(n))

  # Get the geometry type
  geometry_type <- sf::st_geometry_type(shapefile)

  if(!any(geometry_type == "POINT")){
    if(!any(geometry_type == "MULTIPOINT")){
      stop("uc511(getSample) Supplied shapefile must contain POINT or MULTIPOINT geometries.")
    }
  }

  # Check if the geometry type is POINT or MULTIPOINT
  is_point <- base::all(geometry_type == "POINT")
  is_multipoint <- base::all(geometry_type == "MULTIPOINT")

  # if a MULTIPOINT object need to make a POINT object first. should be from $hf.pts.shp.
  if(is_multipoint){
    # get data for the Halton frame.
    base::message("uc511(getSample) is_multipoint.")
    hf_pts <- sf::st_cast(shapefile, "POINT")
    hf_pts <- sf::st_as_sf(hf_pts)
    hf_pts$ID <- base::seq(1, base::length(hf_pts$x))
    hf_pts$uc511SeqID <- seq(1, base::length(hf_pts$x))
    shapefile <- hf_pts
    # now we can make our sample based on randomStart.
  }

  # if is_point TRUE then shapefile from $sample
  if(randomStart){
    # points will already be sorted and each have a valid $uc511SeqID
    base::message("uc511(getSample) is_randomStart.")
    duplicated_pts <- base::rbind(shapefile, shapefile)
    random_start_point <- base::sample(1:length(duplicated_pts$ID), 1)
    sample_indices <- base::seq(random_start_point, (random_start_point + n) - 1, 1)
    sample <- duplicated_pts[sample_indices,]
  }

  # not a random start so just return the first n sample points.
  if(!randomStart){
    # points will already be sorted and each have a valid $uc511SeqID
    # lets sort on $uc511SeqID to be sure...
    base::message("uc511(getSample) is_not_randomStart.")
    shp_sorted <- shapefile[base::order(shapefile$uc511SeqID), ]
    if(n > length(shp_sorted$uc511SeqID)){
      msg <- "uc511(getSample) Warning - n exceeds number of points in shapefile. Returning %s points."
      msgs <- sprintf(msg, length(shp_sorted$uc511SeqID))
      base::message(msgs)
    }
    # get our n sample points.
    sample <- shp_sorted[1:n,]
  }

  # package up objects to be returned.
  result <- base::list(sample = sample)
  return(result)
}


#' @name getHaltonFrame
#'
#' @title Obtain a Halton Frame over a shapefile.
#'
#' @description An internal only function.
#'
#' @details This function was written by Phil Davies.
#'
#' @param shapefile A MULTIPOINT or POINT object that we want to generate a halton frame for.
#' @param J The number of grid cells. A list of 2 values.
#' @param bases Co-prime base for the Halton Sequence.
#' @param i An integer to add to the J parameter elements to expand the Halton Frame in both
#' directions if the required number of sample points cannot be found in the region of interest
#' in the current Halton frame.
#' @param seeds A list of 2 seeds, u1 and u2.
#' @param crs Coordinate reference system for the shapefile.
#'
#' @return A list containing the following variables: hf_, sample, pts.shp, bb.new, seeds
#'
getHaltonFrame <- function(shapefile, J, i, bases, seeds, crs){
  #
  hf_ <- uc511::HaltonFrameBase(J = base::c(J[1]+i, J[2]+i), bases = bases, seeds = seeds)

  # process points returned.
  pts <- hf_$halton_frame
  pts <- base::cbind(base::seq(1, base::dim(pts)[1]), pts)

  # save returned seeds in case they have changed (would only change if initially NULL).
  seeds <- hf_$seeds

  #
  bb <- sf::st_as_sfc(sf::st_bbox(shapefile))
  cntrd <- sf::st_centroid(bb)
  bb.rot <- (bb - cntrd) * uc511::rot(0) + cntrd
  bb.new <- sf::st_as_sfc(sf::st_bbox(bb.rot))

  #
  base::attr(bb.new, "rotation") <- 0
  base::attr(bb.new, "centroid") <- sf::st_coordinates(cntrd)
  pts.shp <- uc511::rotate.scale.coords(coords = pts, bb = bb.new)
  # make sure our shapefile has a CRS (needed for plotting later on).
  sf::st_crs(pts.shp) <- crs
  # always return NULL when just generating a Halton Frame.
  diff_ <- NULL

  result <- base::list(hf_     = hf_,
                       sample  = diff_,
                       pts.shp = pts.shp,
                       bb.new  = bb.new,
                       seeds   = seeds)
  return(result)
}

