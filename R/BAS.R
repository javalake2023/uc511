# BAS.R

#' @import raster
#' @import sf
#' @import Rcpp
#' @useDynLib uc511, .registration = TRUE

#' @name BAS
#'
#' @title Draw spatially balanced samples from areal resources.
#'
#' @description BAS draws spatially balanced samples from areal resources. To draw BAS samples,
#' uc511 requires a study region shapefile and the region’s bounding box. An initial sample size
#' is also needed, which can be easily increased or decreased within uc511 for master sampling
#' applications
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample and later simplified by Phil Davies.
#'
#' @param shapefile Shape file as a polygon (sp or sf) to select sites for.
#' @param n Number of sites to select. If using stratification it is a named vector containing
#' sample sizes of each group.
#' @param boundingbox Bounding box which defines the Master Sample. A bounding box must be
#' supplied.
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
#' @param stratum The name of a column in the data.frame attached to shapefile that defines
#' the strata of interest.
#' @param seeds A vector of 2 seeds, u1 and u2. If not specified, the default is NULL and will
#' be defined randomly using function \code{uc511::generateUVector}.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return A list containing two variables, \code{$sample} containing locations in the BAS sample,
#' in BAS order and \code{$seeds}, the u1 and u2 seeds used to generate the sample.
#'
#' The sample points are returned in the form of a simple feature collection of POINT objects.
#' They have the following attributes:
#' \itemize{
#'   \item \code{SiteID} A unique identifier for every sample point.  This
#'   encodes the BAS order.
#'   \item \code{geometry} The XY co-ordinates of the sample point in the CRS of the original
#'   shapefile.
#'   \item \code{uc511SeqID} A unique identifier for every sample point.  This
#'   encodes the BAS sample order.
#' }

#' @examples
#' \dontrun{
#' data(Fed_MPAs_clipped)
#' # Sample sizes for each stratum:
#' # chose how many sites in each polygon in the dataset
#' N_Zone <- c("Adaptive Management Zone" = 30, "Marine" = 20, "Other" = 40)
#' # Rename the NA value as the function does not accept NA at the moment.
#' Fed_MPAs_clipped$ZONEDESC_E[is.na(Fed_MPAs_clipped$ZONEDESC_E)] <- "Other"
#' # Core Protection is totally within Adaptive Management Zone. Remove it or
#' # make it explicit that it is different.
#' shp.MPAs <- Fed_MPAs_clipped[Fed_MPAs_clipped$ZONEDESC_E != "Core Protection Zone", ]
#' # create a boundingbox
#' bb <- uc511::BoundingBox(shp.MPAs)
#' # Select the Master Sample sites:
#' smp.str <- uc511::BAS(shapefile = shp.MPAs,
#'                       n = N_Zone,
#'                       boundingbox = bb,
#'                       stratum = "ZONEDESC_E",
#'                       verbose = TRUE)
#' plot(sf::st_geometry(shp.MPAs))
#' plot(sf::st_geometry(smp.str$sample), add = T, col= "red", pch = 16)
#' plot(sf::st_geometry(bb), add = T)
#' }
#'
#' @export
BAS <- function(shapefile,
                n = 100,
                boundingbox = NULL,
                panels = NULL,
                panel_overlap = NULL,
                stratum = NULL,
                seeds = NULL,
                verbose = FALSE){

  # validate shapefile and other BAS parameters.
  # validate the shapefile parameter.
  shp_geometry <- sf::st_geometry_type(shapefile)
  if (!base::all(shp_geometry %in% c("MULTIPOLYGON", "POLYGON"))){
    msg <- "uc511(BAS) Unsupported geometry in shapefile, %s."
    msgs <- base::sprintf(msg, shp_geometry)
    base::stop(msgs)
  }

  # A bounding box must be specified.
  if(is.null(boundingbox)){
    msg <- "uc511(BAS) A Bounding Box must be specified. Use uc511::BoundingBox first."
    msgs <- base::sprintf(msg)
    base::stop(msgs)
  }

  # validate panel design if we are using one.
  res <- ValidatePanelDesign(panels, panel_overlap, n)
  panel_design  <- res$panel_design
  number_panels <- res$number_panels
  panel_overlap <- res$panel_overlap
  n             <- res$n

  # stratification wanted?
  if(base::is.null(stratum)){
    # no stratification, just a simple sample wanted.
    result <- uc511::getBASSampleDriver(shapefile = shapefile,
                                        n = n,
                                        bb = boundingbox,
                                        seeds = seeds,
                                        verbose = verbose)
    smp <- result$sample
    seeds <- result$seed
  }else{
    if(base::is.null(base::names(n))) base::stop("uc511(BAS) Need design sample size as n = named vector")
    strata.levels <- base::names(n)

    if(verbose){
      msg <- "uc511(BAS) Stratum: %s.\n"
      msgs <- base::sprintf(msg, strata.levels[1])
      base::message(msgs)
    }
    k.indx <- base::which(shapefile[, stratum, drop = TRUE] == strata.levels[1])
    shp.stratum <- shapefile[k.indx,]
    result <- uc511::getBASSampleDriver(shapefile = shp.stratum,
                                        n = n[1],
                                        bb = boundingbox,
                                        seeds = seeds,
                                        verbose = verbose)
    smp <- result$sample
    seeds <- result$seed
    smp[stratum] <- strata.levels[1]

    if(base::length(n) > 1){
      for(k in 2:base::length(n)){
        if(verbose){
          msg <- "uc511(BAS) Stratum: %s."
          msgs <- base::sprintf(msg, strata.levels[k])
          base::message(msgs)
        }
        k.indx <- base::which(shapefile[, stratum, drop = TRUE] == strata.levels[k])
        shp.stratum <- shapefile[k.indx,]
        result <- uc511::getBASSampleDriver(shapefile = shp.stratum,
                                            n = n[k],
                                            bb = boundingbox,
                                            seeds = seeds,
                                            verbose = verbose)
        smp.s <- result$sample
        seeds <- result$seed
        smp.s[stratum] <- strata.levels[k]
        smp <- base::rbind(smp, smp.s)
      }
    }
  } # end is.null(stratum)

  # go assign panelid's if required.
  res <- uc511::PanelDesignAssignPanelids(smp, panels, panel_overlap, panel_design, number_panels)

  # return the sample and the u1, u2 seeds used.
  result <- base::list(sample = res$sample,
                       seed   = seeds)
  return(result)
}


#' @name getBASSampleDriver
#'
#' @title Title.
#'
#' @description This function repeatedly calls function uc511::getBASSample to generate the BAS
#' sample. Once the requested number of samples within the intersection of the shapefile and the
#' study area have been obtained, the sample and seeds are returned to the caller.
#'
#' @details This function was written by Phil Davies based on origin code by Paul van Dam-Bates
#' from the BASMasterSample package.
#'
#' @param shapefile Shape file as a polygon (sp or sf) to select sites for.
#' @param bb Bounding box which defines the area around the study area. A bounding box must be
#' supplied.
#' @param n Number of sites to select. If using stratification it is a named vector containing
#' sample sizes of each group.
#' @param seeds A vector of 2 seeds, u1 and u2. If not specified, the default is NULL and will
#' be defined randomly.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return A list containing two variables, \code{$sample} containing locations in the BAS sample,
#' in BAS order and \code{$seeds}, the u1 and u2 seeds used to generate the sample.
#'
#' @export
getBASSampleDriver <- function(shapefile, bb, n, seeds, verbose){

  if(is.null(seeds)){
    seeds <- uc511::generateUVector()
  }
  k <- 0

  message("getBASSampleDriver seeds=", seeds)

  # find the first point in the study region (picked at random)
  first.pt <- findFirstStudyRegionPoint(shapefile = shapefile, seeds = seeds)
  # first.pt$first.pt
  # first.pt$seeds
  k <- first.pt$k

  seedshift <- base::c(first.pt$seeds[1] + k - 1, first.pt$seeds[2] + k - 1)

  #pts.sample <- uc511::getBASSample(shapefile = shapefile, bb = bb, n = 2 * n, seeds = seeds)
  #ret_sample <- pts.sample$sample
  #seedshift <- pts.sample$seeds
  #num_samples <- length(ret_sample$SiteID)

  # number of samples required.
  draw <- n * 2
  # just the first point so far, need n.
  num_samples <- 1
  #needed_extra_points <- FALSE
  #
  #seedshift <- first.pt$seeds
  while(num_samples < n){
    #needed_extra_points <- TRUE
    draw <- draw * 2
    last.pt <- num_samples
    endPoint <- last.pt
    pts.sample <- uc511::getBASSample(shapefile = shapefile, bb = bb , n = draw, seeds = seedshift, k = k, endPoint = last.pt)
    #ret_sample <- rbind(ret_sample, pts.sample$sample)
    ret_sample <- pts.sample$sample
    n_samples <- length(ret_sample$SiteID)
    #seedshift <- pts.sample$seeds

    if(verbose){
      msg <- "uc511(getBASSampleDriver) after getBASSample n_samples = %s. k = %s. num_samples = %s\n"
      msgs <- base::sprintf(msg, n_samples, k, endPoint)
      base::message(msgs)
    }
    num_samples <- n_samples #num_samples + n_samples
    #k <- k + 1
  }

  f.pt <- sf::st_cast(first.pt$first.pt, "POINT")
  f.pt <- sf::st_as_sf(f.pt)
  zzz <- sf::st_as_sf(base::data.frame(SiteID = f.pt$ID, f.pt$x))
  ret_sample <- rbind(zzz, ret_sample)
  #sorted_samp <- ret_sample[base::order(ret_sample$SiteID), ]
  sorted_samp <- ret_sample
  sorted_samp$uc511SeqID <- base::seq(1, base::length(sorted_samp$SiteID))

  # return original seeds.
  seeds <- first.pt$seeds

  # return sample and seeds to caller.
  result <- base::list(sample = sorted_samp[1:n,],
                       seeds  = seeds)
  return(result)
}


#' @export
getBASSample <- function(shapefile, bb, n, seeds, k = 0, endPoint = 0){

  first_point_in_region <- FALSE
  seedshift <- seeds

  # Scale and shift Halton to fit into bounding box
  bb.bounds <- sf::st_bbox(bb)
  scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
  shift.bas <- bb.bounds[1:2]

  # count number times we have to try and find first BAS point in the study region.
  attempts_to_find_first_point <- 0

  #while(!first_point_in_region){
    #pts <- uc511::cppRSHalton(n = draw, seeds = seedshift, bases = c(2, 3), boxes = halt.rep, J = J)
    #pts <- pts[1:draw,]
    res <- uc511::cppRSHalton_br(n = n, seeds = seedshift, bases = base::c(2, 3))
    pts <- res$pts
    siteid <- base::seq(from = 1, to = n, by = 1)
    pts <- base::cbind(siteid, pts)

    xy <- base::cbind(pts[,2]*scale.bas[1] + shift.bas[1], pts[,3]*scale.bas[2] + shift.bas[2])

    pts.coord <- sf::st_as_sf(base::data.frame(SiteID = pts[,1], xy), coords = c(2, 3))

    sf::st_crs(pts.coord) <- sf::st_crs(bb)
    # find the intersection.
    pts.intersect <- pts.coord[shapefile,]
    #if(length(pts.intersect$SiteID) == 0){
      # no points in intersection, so keep trying...
    #  attempts_to_find_first_point <- attempts_to_find_first_point + 1
    #  seedshift <- uc511::generateUVector()
    #  next
    #}
    # is the first point of pts.coord in the shapefile?
    #if(pts.coord$SiteID[1] == pts.intersect$SiteID[1]){
    #  first_point_in_region <- TRUE
    #} else{
    #  attempts_to_find_first_point <- attempts_to_find_first_point + 1
    #  seedshift <- uc511::generateUVector()
    #}
  #} # end while

  #if(attempts_to_find_first_point > 0){
  #  msg <- "uc511() First point in study area found in %s attempts."
  #  msgs <- base::sprintf(msg, attempts_to_find_first_point)
  #  base::message(msgs)
  #}

  result <- base::list(sample = pts.intersect,
                       seeds  = seedshift)
  return(result)
}


#' @export
generateUVector <- function(){
  u1 <- base::floor(stats::runif(1, 0, 2^11))
  u2 <- base::floor(stats::runif(1, 0, 3^7))
  seeds <- c(u1, u2)
  return(seeds)
}


# 1. Set J1 = 4 and J2 = 3.
# 2. Generate B = 2^J1 x 3^J2 points from a random-start Halton sequence H
#    with a random seed (u1, u2).
# 3. Find points from H in the study area. Call this set S. If S is empty,
#    increment J1 and J2 and go to step 2.
# 4. Randomly choose a point from S. Let xk be this point where k is the 'site index'
#    (I think that's what we call it).
# 5. Set the seed to (u1 + k - 1, u2 + k -1).

# For example, let (u1, u2) = (1, 5) and S = {x2, x6, x7}.
# If x6 is randomly chosen, then the new seed is (1 + 6 - 1, 5 + 6 - 1) = (6, 10)
# (the sixth point in H).

# The only difference is that the random-start Halton sequence must be length B.

#' @export
findFirstStudyRegionPoint <- function(shapefile, seeds){
  n <- (2^4) * (3^3)

  if(is.null(seeds)){
    seeds <- uc511::generateUVector()
  }

  J <- base::c(4, 3)
  i <- 0
  bases <- base::c(2, 3)
  crs <- sf::st_crs(shapefile)

  result <- getHaltonFrame(shapefile, J, i, bases, seeds, crs)
  hf_ <- result$hf_
  diff_ <- result$sample     # will always be NULL here.
  pts.shp <- result$pts.shp
  bb.new <- result$bb.new
  seeds <- result$seeds

  pts <- hf_$halton_frame
  tmp <- sf::st_cast(pts.shp, "POINT")
  tmp <- sf::st_as_sf(tmp)
  tmp$ID <- base::seq(1, base::dim(pts)[1])
  diff_ <- sf::st_intersection(tmp, shapefile)

  # find number of points within our shapefile.
  pts_in_intersection <- base::length(sf::st_cast(sf::st_union(diff_), "POINT"))

  # will need code here to detect 0 points in intersection and redrive.

  msg <- "uc511(findFirstStudyRegionPoint) Points in intersection: %s."
  msgs <- sprintf(msg, pts_in_intersection)
  base::message(msgs)

  s <- diff_

  base::set.seed(seeds[1] + seeds[2])
  k <- base::sample(pts_in_intersection, 1)
  message("findFirstStudyRegionPoint k=", k)

  first.pt <- s[k,]
  # Set the seed to (u1 + k - 1, u2 + k -1)
  #seeds <- base::c(seeds[1] + k - 1, seeds[2] + k - 1)

  result <- base::list(first.pt = first.pt,
                       seeds    = seeds,       # original seeds here!
                       k        = k)
  return(result)
}

