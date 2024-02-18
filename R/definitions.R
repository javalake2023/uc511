#' @name BoundingBox
#'
#' @title Build a new Master Sample with a random rotation and seed.
#'
#' @description Randomly generate a seed from 10,000 possible values in right now 2 dimensions.
#' Note that in van Dam-Bates et al. (2018) we required that the random seed falls into main object shape, such
#' as one of the islands in New Zealand, or within marine environment for BC west coast. However, with a random rotation,
#' we are able to ignore that detail. If this function is used without a random rotation, we recommend running it until
#' the first master sample point does indeed fall within the largest scale of the master sample use.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample and later ported to this package, uc511.
#'
#' @param shp Spatial feature that defines the boundary of the area to define a bounding box over.
#' @param d Dimension of the new Master Sample, at this stage we only work with d=2.
#' @param showOutput Print the rotation and random seed when it is generated.
#' @param rotate Boolean of whether or not to randomly rotate the bounding box.
#'
#' @return bounding box for a master sample.
#'
#' @examples
#' \dontrun{
#' data(NS_bioregion)
#' bb <- uc511::BoundingBox(shp = NS_bioregion)         # Vertically aligned master sample bounding box.
#' bb.rot <- uc511::rotate.shp(bb, bb, back = TRUE) # Actual bounding box.
#' plot(sf::st_geometry(NS_bioregion))
#' plot(sf::st_geometry(bb.rot), add = TRUE)
#' }
#'
#' @export
BoundingBox <- function(shp, d = 2, showOutput = TRUE, rotate = FALSE)
{
  # generate 2 seed values.
  seed <- base::floor(stats::runif(d, 0, 10000))
  # We always use base 2,3
  base <- c(2, 3, 5)[1:d]

  if(rotate) {
    theta <- stats::runif(1, -base::pi, base::pi)
  }else{
    theta <- 0
  }

  # Just work with sf:
  if (class(shp)[1] != "sf"){
    shp <- sf::st_as_sf(shp)
  }

  # Create a Random Rotation:
  build.bb <- rotate.bb(shp, theta = theta)
  sf::st_crs(build.bb) <- sf::st_crs(shp)
  attr(build.bb, "seed") <- seed

  if(showOutput){
    #msg <- "uc511(BoundingBox) Seed: %s.\n"
    #msgs <- sprintf(msg, seed)
    #base::message(msgs)

    cat("uc511(BoundingBox) Seed:", seed, ".", "\n")

    msg <- "uc511(BoundingBox) Rotation: %s Radians.\n"
    msgs <- sprintf(msg, theta)
    base::message(msgs)
  }
  return(build.bb)
}


#' @name getProj
#'
#' @title Define spatial objects in projection of the master sample.
#'
#' @description Default projection of the master sample. Needed for consistency for the entire bounding box.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample and later ported to this package, uc511.
#'
#' @return Spatial object projection.
#'
#' @export
getProj <- function()
{	#BC Albers
  #http://spatialreference.org/ref/epsg/3005/
  msproj <- "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  return(msproj)
}


#' @name getBB
#'
#' @title Get the bounding box for other functions.
#'
#' @description Bounding box under NAD83/BC Albers.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample and later ported to this package, uc511.
#'
#' @return The bounding box.
#'
#' @export
getBB <- function()
{
  bb.df <- c("xmin" = 85148, "ymin" = 33745, "xmax" = 1280999, "ymax" = 1351981)
  bb <- sf::st_as_sfc(sf::st_bbox(bb.df))

  attr(bb, "rotation") <- 0
  attr(bb, "centroid") <- sf::st_centroid(bb)

  sf::st_crs(bb) <- sf::st_crs(getProj())
  return(bb)
}


#' @name getSeed
#'
#' @title Random seed definition for Western Canada Marine Master Sample.
#'
#' @description Defines the random seed specific to the Western Canada Marine Master Sample.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample and later ported to this package, uc511.
#'
#' @return The random seeds for the Western Canada Marine Master Sample.
#'
#' @export
getSeed <- function()
{
  seed <- c(37916, 85846)
  return(seed)
}


#' @name getRotation
#'
#' @title Returns the radians of rotation of the bounding box.
#'
#' @description Defines the random rotation specific to the Western Canada Marine Master Sample.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample and later ported to this package, uc511.
#'
#' @return The radians of rotation of the bounding box.
#'
#' @export
getRotation <- function()
{
  return(0)
}
