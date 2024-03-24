# BAS.R

#' @name BAS
#'
#' @title Select points from a polygon using a BAS Master Sample.
#'
#' @description This is the main function for selecting sites using the BAS master
#' sample. It assumes that you have already defined the master sample using the
#' BoundingBox() function.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
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
#' be defined randomly.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return A list containing two variables, \code{$sample} containing locations in the BAS sample,
#' in BAS order and \code{$seeds}, the u1 and u2 seeds used to generate the sample.
#'
#' The sample points are returned in the form of a simple feature collection of POINT objects.
#' They have the following attributes:
#' \itemize{
#'   \item \code{SiteID}: A unique identifier for every sample point.  This
#'   encodes the BAS order.
#'   \item \code{geometry}: The XY co-ordinates of the sample point in the CRS of the original
#'   shapefile.
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
    # no stratification, just a simple sample.
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
#' @title Select points from a polygon using a BAS Master Sample.
#'
#' @description This function repeatedly calls function uc511::getBASSample to generate the BAS
#' sample. Once the requested number of samples within the intersection of the shapefile and the
#' study area have been obtained, the sample and seeds are returned to the caller.
#'
#' @details This function was written by Phil Davies based on origin code by Paul van Dam-Bates
#' from the BASMasterSample package.
#'
#' @param shapefile Shape file as a polygon (sp or sf) to select sites for.
#' @param n Number of sites to select. If using stratification it is a named vector containing
#' sample sizes of each group.
#' @param boundingbox Bounding box which defines the Master Sample. A bounding box must be
#' supplied.
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
    seeds <- c(0, 0)
  }
  k <- 0

  pts.sample <- uc511::getBASSample(shapefile = shapefile, bb = bb, n = n, seeds = seeds)
  ret_sample <- pts.sample$sample
  num_samples <- length(ret_sample$SiteID)
  # number of samples required.
  draw <- n
  seedshift <- seeds
  while(num_samples < n){
    draw <- draw * 2
    last.pt <- num_samples
    endPoint <- last.pt
    pts.sample <- uc511::getBASSample(shapefile = shapefile, bb = bb , n = draw, seeds = seedshift, k = k, endPoint = last.pt)
    ret_sample <- pts.sample$sample
    n_samples <- length(ret_sample$SiteID)
    seedshift <- pts.sample$seeds

    if(verbose){
      msg <- "uc511(getBASSampleDriver) after getSample %s. k = %s. endPoint = %s\n"
      msgs <- base::sprintf(msg, num_samples, k, endPoint)
      base::message(msgs)
    }
    num_samples <- n_samples
    k <- k + 1
  }
  # return sample and seeds to caller.
  result <- base::list(sample = pts.sample$sample[1:n,],
                       seeds  = seeds)
  return(result)
}


#' @export
getBASSample <- function(shapefile, bb, n, seeds, k = 0, endPoint = 0){

  seedshift <- seeds

  # Scale and shift Halton to fit into bounding box
  bb.bounds <- sf::st_bbox(bb)
  scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
  shift.bas <- bb.bounds[1:2]

  #pts <- uc511::cppRSHalton(n = draw, seeds = seedshift, bases = c(2, 3), boxes = halt.rep, J = J)
  #pts <- pts[1:draw,]
  res <- uc511::cppRSHalton_br(n = n, seeds = seedshift, bases = base::c(2, 3))
  pts <- res$pts
  siteid <- base::seq(from = 1, to = n, by = 1)
  pts <- base::cbind(siteid, pts)

  xy <- base::cbind(pts[,2]*scale.bas[1] + shift.bas[1], pts[,3]*scale.bas[2] + shift.bas[2])

  pts.coord <- sf::st_as_sf(base::data.frame(SiteID = pts[,1] + endPoint, xy), coords = c(2, 3))

  sf::st_crs(pts.coord) <- sf::st_crs(bb)
  # find the intersection.
  pts.coord <- pts.coord[shapefile,]

  result <- base::list(sample = pts.coord,
                       seeds  = seedshift)
  return(result)
}

