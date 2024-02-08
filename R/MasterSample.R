########################
#Master Sample Code
#Paul van Dam-Bates
########################

#' @import raster
#' @import sf
#' @import Rcpp
#' @useDynLib uc511, .registration = TRUE


#' @name masterSampleSelect
#'
#' @title Generate sample points using BAS master sample.
#'
#' @description This function has no features and is generally expected to be run by the wrapper function.
#' masterSample that ensures all the correct pieces are passed and adds the additional feature of stratification.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
#'
#' @param shp Shape file as a polygon to select sites for.
#' @param n Number of sites to select.
#' @param bb Bounding box which defines the Master Sample. Default is BC Marine Master Sample.
#' @param nExtra An efficiency problem of how many extra samples to draw before spatial clipping to shp.
#' @param printJ Boolean if you want to see J, how many cuts of space are required to generate the sample efficiently.
#' @param inclSeed Boolean need something here
#'
#' @returns Sample points.
#'
#' @export
masterSampleSelect <- function(shp, n = 100, bb = NULL, nExtra = 0, printJ = FALSE, inclSeed = NULL){

  # future option: optInclProb - return include probabilities or not.
  optInclProb <- FALSE

  # We always use base 2,3
  base <- c(2,3)

  if(is.null(inclSeed)) inclSeed <- base::floor(stats::runif(1,1,10000))

  # Updating to work for sf only. Start here...
  if (class(shp)[1] != "sf"){
    shp <- sf::st_as_sf(shp)
  }

  # Set up Western Canada Marine Master Sample as default, general for any.
  # bb now includes its rotation as well.
  if(is.null(bb)){
    bb <- getBB()
    msproj <- getProj()
    seed <- getSeed()
  }else{
    msproj <- sf::st_crs(bb)$proj4string
    seed <- base::attr(bb, "seed")[1:2]	# 3rd dimension is not yet supported...
  }

  orig.crs <- NULL
  if(sf::st_crs(shp) != sf::st_crs(msproj)){
    orig.crs <- sf::st_crs(shp)$proj4string
    shp <- sf::st_transform(shp, msproj)
  }

  #Scale and shift Halton to fit into bounding box
  bb.bounds <- sf::st_bbox(bb)
  scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
  shift.bas <- bb.bounds[1:2]

  cntrd <- base::attr(bb, "centroid")
  theta <- base::attr(bb, "rotation")

  # We can use Halton Boxes to speed up code when the polygons are small and all over the place.
  # Kind of like magic!
  draw <- n + nExtra

  J <- c(0, 0)
  # Tricky here to figure out where the shape is on the vertical Halton box.
  shp.rot <- rotate.shp(shp, bb, back = FALSE)
  hal.frame <- shape2Frame(shp.rot, J = J, bb = bb, projstring = msproj)
  area.shp <- as.numeric(base::sum(sf::st_area(shp)))
  # Subset again:
  while(area.shp < 0.25*as.numeric(sf::st_area(hal.frame))[1]){
    if(base[2]^J[2] > base[1]^J[1]){
      J[1] <- J[1] + 1
    }else{
      J[2] <- J[2] + 1
    }
    hal.frame <- shape2Frame(shp.rot, J = J, bb = bb, projstring = msproj)
  }

  hal.fr2 <- uc511::rotate.shp(hal.frame, bb)
  hal.indx <- which(raster::rowSums(sf::st_intersects(hal.fr2, shp, sparse = FALSE)) > 0)
  # need to remove prior to submitting to CRAN.
  if(length(hal.indx) < 2){browser()}

  #hal.pts <- (sf::st_centroid(hal.frame) %>% sf::st_coordinates())[hal.indx,]	## Hack to react to changes in sf. Need to get coordinates and then subset now. Fix it better in the future.
  hal.tmp <- sf::st_centroid(hal.fr2) #hal.frame ???? #[hal.indx,]
  hal.pts <- sf::st_coordinates(hal.tmp)[hal.indx,]

  # Find the corner Halton Pts
  box.lower <- t(base::apply(hal.pts, 1, FUN = function(x){(x - shift.bas)/scale.bas}))
  A <- uc511::GetBoxIndices(box.lower, base, J)
  halt.rep <- uc511::SolveCongruence(A, base, J)
  B <- base::prod(c(2, 3) ^ J)

  # We like to know how many divisions we had to make...
  if(printJ){
    msg <- "uc511(masterSampleSelect) Number of divisions made (J=) %s.\n"
    msgs <- base::sprintf(msg, J)
    base::message(msgs)
  }

  getSample <- function(k = 0, endPoint = 0){
    print("getSample")
    seed <- c(seed, inclSeed)
    if(k == 0){ seedshift <- seed
    }else {
      seedshift <- endPoint + seed
    }
    print("seedshift:")
    print(seedshift)
    #pts <- uc511::cppRSHalton(n = draw, seeds = seedshift, bases = c(2, 3, 5), boxes = halt.rep, J = J)
    #pts <- pts[1:draw,]
    res <- uc511::cppRSHalton_br(n = draw, seeds = seedshift, bases = c(2, 3, 5))
    pts <- res$pts
    siteid <- base::seq(from = 1, to = draw, by = 1)
    pts <- base::cbind(siteid, pts)

    xy <- base::cbind(pts[,2]*scale.bas[1] + shift.bas[1], pts[,3]*scale.bas[2] + shift.bas[2])
    if(theta != 0) xy <- sweep ( sweep(xy, 2,  cntrd, FUN = "-") %*% uc511::rot(-theta), 2,  cntrd, FUN = "+")
    if(optInclProb){
      pts.coord <- sf::st_as_sf(data.frame(SiteID = pts[,1] + endPoint, xy, inclProb = pts[,4]), coords = c(2,3))
    } else {
      pts.coord <- sf::st_as_sf(data.frame(SiteID = pts[,1] + endPoint, xy), coords = c(2,3))
    }
    sf::st_crs(pts.coord) <- sf::st_crs(bb)
    pts.coord <- pts.coord[shp,]
    return(pts.coord)
  }

  pts.sample <- getSample()
  while(nrow(pts.sample) == 0){
    draw <- draw * 2
    pts.sample <- getSample()
  }

  di <- 1
  while(nrow(pts.sample) < n){
    last.pt <- pts.sample$SiteID[base::nrow(pts.sample)]
    new.pts <- getSample(k = di, endPoint = last.pt)
    if(nrow(new.pts) > 0) pts.sample <- base::rbind(pts.sample, new.pts)
    di <- di + 1
  }

  smp <- pts.sample[1:n,]
  if(!is.null(orig.crs)){
    smp <- sf::st_transform(smp, orig.crs)
  }
  return(smp)
}


#' @name getBASMasterSample
#'
#' @title Select points from a polygon using a BAS Master Sample.
#'
#' @description This is the main function for selecting sites using the BAS master
#' sample. It assumes that you have already defined the master sample using the
#' buildMS() function or will be selecting a marine master sample site in BC.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
#'
#' @param shp Shape file as a polygon (sp or sf) to select sites for.
#' @param n Number of sites to select. If using stratification it is a named vector containing sample sizes of each group.
#' @param bb Bounding box which defines the Master Sample. Default is BC Marine Master Sample.
#' @param stratum Name of column of data.frame attached to shapefile that makes up the strata.
#' @param nExtra An efficiency problem of how many extra samples to draw before spatial clipping to shp.
#' @param quiet Boolean if you want to see any output printed to screen. Helpful if taking a long time.
#' @param inclSeed Boolean need something here
#'
#' @return A master sample.

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
#' # Select the Master Sample sites:
#' smp.str <- uc511::getBASMasterSample(shp.MPAs, N = N_Zone, stratum = "ZONEDESC_E", quiet = FALSE)
#' plot(sf::st_geometry(shp.MPAs))
#' plot(sf::st_geometry(smp.str), add = T, col= "red", pch = 16)
#' }
#'
#' @export
getBASMasterSample <- function(shp, n = 100, bb = NULL, stratum = NULL, nExtra = 0, quiet = FALSE, inclSeed = NULL)
{
  if(is.null(inclSeed)) inclSeed <- base::floor(stats::runif(1,1,10000))
  if(is.null(stratum)){
    smp <- masterSampleSelect(shp = shp, n = n, bb = bb, nExtra = nExtra, inclSeed = inclSeed)
  }else{
    if(is.null(names(n))) return("Need design sample size as n = named vector")
    strata.levels <- names(n)

    if(!quiet){
      msg <- "uc511(getBASMasterSample) Stratum: %s.\n"
      msgs <- base::sprintf(msg, strata.levels[1])
      base::message(msgs)
    }
    k.indx <- which(shp[, stratum, drop = TRUE] == strata.levels[1])
    shp.stratum <- shp[k.indx,] #%>% st_union()	# ? Not sure if this is necessary... slowed things down too much!
    smp <- masterSampleSelect(shp.stratum, n = n[1], bb = bb, nExtra = nExtra, printJ = !quiet, inclSeed)
    smp[stratum] <- strata.levels[1]

    if(base::length(n) > 1){
      for(k in 2:base::length(n)){
        if(!quiet){
          msg <- "uc511(getBASMasterSample) Stratum: %s."
          msgs <- base::sprintf(msg, strata.levels[k])
          base::message(msgs)
        }
        k.indx <- which(shp[, stratum, drop = TRUE] == strata.levels[k])
        shp.stratum <- shp[k.indx,] ## %>% st_union()	# Needed?
        smp.s <- masterSampleSelect(shp = shp.stratum, n = n[k], bb = bb, nExtra = nExtra, printJ = !quiet, inclSeed = inclSeed)
        smp.s[stratum] <- strata.levels[k]
        smp <- base::rbind(smp, smp.s)
      }
    }
  } # end is.null(stratum)
  result <- base::list(sample = smp,
                       seed   = inclSeed)
  return(result)
}


#' @name point2Frame
#'
#' @title Get the Halton Boxes around a point resource ordered by BAS Master Sample.
#'
#' @description Get the Halton Boxes around a point resource ordered by BAS Master Sample.
#' If you pass discrete points to this function it will return how the Halton
#' frame is represented and what the ordering is.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
#'
#' @param pts Shape file as a polygon (sp or sf) to generate the boxes for.
#' @param bb Bounding box which defines the Master Sample. Default is BC Marine Master Sample.
#' @param base Co-prime base of Halton sequence
#' @param J Definition for the number of grid cells of Halton frame.
#' @param size Physical target size of Halton boxes (square ish) if J is NULL.
#'
#' @return how the Halton frame is represented and what the ordering is.

#' @examples
#' \dontrun{
#' # If you haven't already installed this:
#' install.packages('bcmapsdata', repos='https://bcgov.github.io/drat/')
#' library(bcmapsdata)
#' library(bcmaps)
#' cities <- bcmaps::get_layer("bc_cities")
#' bb <- uc511::buildMS(cities, d = 2, FALSE)
#' # For visibility will make boxes 10 km
#' cities.halton <- uc511::point2Frame(pts = cities, bb = bb, size = 10000)
#' plot(sf::st_geometry(cities), pch = 20, cex = 0.1)
#' plot(sf::st_geometry(cities.halton), add= TRUE)
#' # What is the actual area of a Halton box?
#' sf::st_area(cities.halton[1,])/1000^2
#' }
#'
#' @export
point2Frame <- function(pts, bb = NULL, base = c(2,3), J = NULL, size = 100)
{
  # Updating to work for sf only. Start here...
  if (class(pts)[1] != "sf"){
    pts <-  sf::st_as_sf(pts)
  }

  # Set up Western Canada Marine Master Sample as default, general for any.
  # bb now includes its rotation as well.
  if(is.null(bb)){
    bb <- getBB()
    msproj <- getProj()
    seed <- getSeed()
  }else{
    msproj <- sf::st_crs(bb)$proj4string
    seed <- base::attr(bb, "seed")[1:2]	# 3rd dimension is not yet supported...
  }

  orig.crs <- NULL
  if(sf::st_crs(pts) != sf::st_crs(msproj)){
    orig.crs <- sf::st_crs(pts)$proj4string
    pts <- sf::st_transform(pts, msproj)
  }

  #Scale and shift Halton to fit into bounding box
  bb.bounds <- sf::st_bbox(bb)
  scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
  shift.bas <- bb.bounds[1:2]

  cntrd <- base::attr(bb, "centroid")
  theta <- base::attr(bb, "rotation")

  if(is.null(J)) J <- base::ceiling(base::log(scale.bas/size)/base::log(base))

  xy.r <- rotate.shp(pts, bb, back = FALSE)
  xy <- sf::st_coordinates(xy.r)
  xy <- base::cbind((xy[,1] - shift.bas[1])/scale.bas[1], (xy[,2] - shift.bas[2])/scale.bas[2])
  lx <- base::floor(xy[,1] / (1/base[1]^J[1]))/(base[1]^J[1])
  ly <- base::floor(xy[,2] / (1/base[2]^J[2]))/(base[2]^J[2])
  A <- GetBoxIndices(base::cbind(lx,ly), base, J)

  #Only build the boxes for unique ones.
  frame.order <- uc511::SolveCongruence(A, base, J)

  dupli <- base::duplicated(frame.order)

  lx <- lx*scale.bas[1] + shift.bas[1]
  ly <- ly*scale.bas[2] + shift.bas[2]

  # Adjust everything for the Master Sample random seed.
  a1 <- seed[1:2] %% base^J
  boxInit <- uc511::SolveCongruence(base::matrix(a1, ncol = 2), base = base, J = J)
  B <- base::prod(base^J)
  # Adjusted index:
  frame.order <- ifelse( frame.order < boxInit, B + ( frame.order - boxInit),  frame.order - boxInit)

  polys <- base::cbind("xmin" = lx, "ymin" = ly, "xmax" = lx + scale.bas[1]/base[1]^J[1], "ymax" = ly + scale.bas[2]/base[2]^J[2])
  halt.boxes <- do.call("c", apply(polys[!dupli,] , 1, FUN = function(x){sf::st_as_sfc(sf::st_bbox(x))}))
  sf::st_crs(halt.boxes) <- sf::st_crs(bb)
  halt.rot <- rotate.shp(halt.boxes, bb, back = TRUE)
  shp.out <- sf::st_sf(halt.rot, data.frame(HaltonIndex = frame.order[!dupli]))

  return(shp.out)
}


#' @name getIndividualBoxIndices
#'
#' @title Get the Halton Box Index for a point.
#'
#' @description If you pass discrete points to this function it will return the Halton Index according to the
#' Master Sample as defined by the bounding box. If the points are closer together than the defined box
#' size then they are treated as the same sample unit. We do not expect this to be spatially balanced and do
#' not recommend using this function alone to select a sample.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
#'
#' @param pts Shape file as a polygon (sp or sf) to generate the boxes for.
#' @param J Definition for the number of grid cells of Halton frame.
#' @param bb Bounding box which defines the Master Sample. Default is BC Marine Master Sample.
#' @param size Physical target size of Halton boxes (square ish) if J is NULL.
#'
#' @return Return the Halton Index for all "pts" that are passed to this function.

#' @examples
#' \dontrun{
#' library(bcmaps)
#' # If you haven't already installed this:
#' install.packages('bcmapsdata', repos='https://bcgov.github.io/drat/')
#' library(bcmapsdata)
#' cities <- bcmaps::get_layer("bc_cities")
#' bb <- uc511::buildMS(cities, d = 2, FALSE)
#' # For visibility will make boxes 10 km
#' cities.ord <- uc511::getIndividualBoxIndices(pts = cities, bb = bb, size = 100)
#' plot(st_geometry(cities), pch = 20)
#' plot(st_geometry(cities.ord[rank(cities.ord$HaltonIndex) < 15,]), add= TRUE, col = "red", cex = 1.5)
#' }
#'
#' @export
getIndividualBoxIndices <- function(pts, J = NULL, bb, size = 100)
{
  # We always use base 2,3
  base <- c(2,3)

  # Updating to work for sf only. Start here...
  if (class(pts)[1] != "sf"){
    pts <-  sf::st_as_sf(pts)
  }

  # Set up Western Canada Marine Master Sample as default, general for any.
  # bb now includes its rotation as well.
  if(is.null(bb)){
    bb <- getBB()
    msproj <- getProj()
    seed <- getSeed()
  }else{
    msproj <- sf::st_crs(bb)$proj4string
    seed <- base::attr(bb, "seed")[1:2]	# 3rd dimension is not yet supported...
  }

  orig.crs <- NULL
  if(sf::st_crs(pts) != sf::st_crs(msproj)){
    orig.crs <- sf::st_crs(pts)$proj4string
    pts <- sf::st_transform(pts, msproj)
  }

  # Scale and shift Halton to fit into bounding box
  bb.bounds <- sf::st_bbox(bb)
  scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
  shift.bas <- bb.bounds[1:2]

  cntrd <- base::attr(bb, "centroid")
  theta <- base::attr(bb, "rotation")

  if(is.null(J)) J <- base::ceiling(base::log(scale.bas/size)/base::log(base))

  B <- base::prod(base^J)

  Bx <- base[1]^J[1]
  By <- base[2]^J[2]

  xy <- sf::st_coordinates(pts)
  # Rotate pts to the bounding box:
  if(theta != 0) xy <- base::sweep(base::sweep(xy, 2,  cntrd, FUN = "-") %*% rot(theta), 2,  cntrd, FUN = "+")
  # Scale to 0-1
  xy <- base::cbind((xy[,1] - shift.bas[1])/scale.bas[1], (xy[,2] - shift.bas[2])/scale.bas[2])
  Axy <- base::cbind(floor((xy[,1] + 2*.Machine$double.eps)*Bx), base::floor((xy[,2] + 2*.Machine$double.eps)*By))

  haltonIndex <- uc511::SolveCongruence(Axy, base = base, J = J)

  # Adjust everything for the Master Sample random seed.
  a1 <- seed[1:2] %% base^J
  boxInit <- uc511::SolveCongruence(base::matrix(a1, ncol = 2), base = base, J = J)

  # Adjusted index:
  haltonIndex <- ifelse(haltonIndex < boxInit, B + (haltonIndex - boxInit), haltonIndex - boxInit)
  # Return the Halton Index for all "pts" in dat that are passed to this function.

  pts$HaltonIndex <- haltonIndex
  return(pts)
}
