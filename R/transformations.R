# transformations.R

#' @name rot
#'
#' @title Generate a rotation matrix for rotating objects later.
#'
#' @description Generate a rotation matrix for rotating objects later.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
#'
#' @param a radians of rotation.
#'
#' @return Matrix
#'
#' @examples
#' \dontrun{
#' rot(pi)
#' }
#'
#' @export
rot <- function(a){
  base::matrix(c(base::cos(a), base::sin(a), -base::sin(a), base::cos(a)), 2, 2)
}


#' @name rotate.bb
#'
#' @title Rotate Bounding box by theta radians
#'
#' @description Given some shp defined as the boundary of interest, rotate it around the centroid
#' and return the rotation and the centroid as attributes. This is used for
#' defining a Master Sample bounding box that has random rotation while ensuring that
#' the new rotated bounding box fits the shp.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
#'
#' @param shp A spatial file with the spatial boundary of the sample.
#' @param theta Radians of rotation. Positive to the right of pi/2, negative to the left.
#'
#' @return Matrix
#'
#' @examples
#' \dontrun{
#' data(NS_bioregion)
#' bb.new <- uc511::rotate.bb(NS_bioregion, -pi/3)
#' }
#'
#' @export
rotate.bb <- function(shp, theta){
  if (base::class(shp)[1] != "sf"){
    shp <- sf::st_as_sf(shp)
  }

  bb <- sf::st_as_sfc(sf::st_bbox(shp))

  cntrd <- sf::st_centroid(bb)
  bb.rot <- (bb - cntrd) * uc511::rot(theta) + cntrd
  bb.new <- sf::st_as_sfc(sf::st_bbox(bb.rot))

  base::attr(bb.new, "rotation") = theta
  base::attr(bb.new, "centroid") = sf::st_coordinates(cntrd)
  return(bb.new)
}


#' @name rotate.scale.coords
#'
#' @title Scale and rotate points from the unit square to a defined projection.
#'
#' @description Given some coordinates on \[0,1)x\[0,1), shift and scale them to the bounding box, and then rotate
#' them given the bounding box rotation defined by the Master Sample.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
#'
#' @param coords Output from RSHalton() to be converted to the spatial surface of interest.
#' @param bb Special shape file defining the bounding box with attributes for centroid and rotation.
#' @param back Boolean for whether or not the rotation is back to the original rotated bounding box.
#'
#' @return sf spatial points with projection defined in bb.
#'
#' @examples
#' \dontrun{
#' pts <- uc511::cppRSHalton(n = 10)
#' bb <- uc511::getBB()
#' pts.shp <- uc511::rotate.scale.coords(coords = pts, bb)
#' }
#'
#' @export
rotate.scale.coords <- function(coords, bb, back = TRUE){

  coords <- coords[, 2:3]
  theta <- ifelse(back, -1, 1) * base::attr(bb, "rotation")	# Rotate backwards
  cntrd <- base::attr(bb, "centroid")

  bb.bounds <- sf::st_bbox(bb)
  bb.scale <- base::diag(2) * (bb.bounds[3:4] - bb.bounds[1:2])

  coords.tmp <- sf::st_multipoint(coords, dim = "XY")
  coords <- sf::st_geometry(coords.tmp)
  coords.scale <- coords*bb.scale + bb.bounds[1:2]
  coords.rot <- (coords.scale - cntrd) * uc511::rot(theta) + cntrd

  sf::st_crs(coords.rot) <- sf::st_crs(bb)
  return(coords.rot)
}


#' @name rotate.shp
#'
#' @title Rotate a polygon around the centroid of a Master Sample bounding box.
#'
#' @description Given some polygon within the bounding box of a Master Sample rotate it by theta defined
#' by that bounding box either backwards or forwards.
#'
#' @param shp Any polygon within the sample frame defined as a spatial features object.
#' @param bb Special shape file defining the bounding box with attributes for centroid and rotation.
#' @param back Boolean for whether or not the rotation is back to the original rotated bounding box.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
#'
#' @return rotated sf spatial object.
#'
#' @examples
#' \dontrun{
#' data(NS_bioregion)
#' bb <- uc511::getBB()
#' pts.shp <- uc511::rotate.shp(shp = NS_bioregion, bb = bb)
#' }
#'
#' @export
rotate.shp <- function(shp, bb, back = TRUE){
  theta <- ifelse(back, -1, 1) * base::attr(bb, "rotation")	# Rotate backwards
  cntrd <- base::attr(bb, "centroid")

  shp <- sf::st_transform(shp, sf::st_crs(bb))
  shp <- sf::st_geometry(shp)
  shp.rot <- (shp - cntrd) * uc511::rot(theta) + cntrd
  sf::st_crs(shp.rot) <- sf::st_crs(bb)
  return(shp.rot)
}
