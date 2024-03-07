# HaltonFunctions.R

#' @name makeFrame
#'
#' @title Create a Halton Frame raster, defined in sf based on the bounding box master sample.
#'
#' @description Make a Halton Frame based on B = 2^J\[1\]*3^J\[2\] grid cells. If rotation is required, will return rotated.
#' This function is an internal function simply to select sub BAS points without having to do spatial clipping at the point level.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
#'
#' @param base Co-prime base for BAS, do not change from 2,3.
#' @param J Definition for the number of grid cells of Halton frame.
#' @param bb Bounding box shapefile with centroid, random seed, rotation.
#' @param rotate Boolean if you want to rotate shape before exporting.
#'
#' @return rotated sf spatial object.
#'
#' @examples
#' \dontrun{
#' bb <- uc511::getBB()
#' haltonFrame <- uc511::makeFrame(J = c(8,4), bb = bb)
#' }
#'
#' @export
makeFrame <- function(base = c(2,3), J = c(2,2), bb, rotate = FALSE){

  b.bounds <- sf::st_bbox(bb)
  B <- base::prod(base^J)
  halt.grid <- raster::raster(raster::extent(base::matrix( b.bounds, 2, 2 )), nrow=base[2]^J[2], ncol=base[1]^J[1])
  halt.grid <- raster::rasterToPolygons(halt.grid)
  raster::projection(halt.grid) <- sf::st_crs(bb)$proj4string
  if(rotate) return(uc511::rotate.shp(sf::st_as_sf(halt.grid), bb))
  return(sf::st_as_sf(halt.grid))
}


#' @name shape2Frame
#'
#' @title Clip a Halton Frame based on the bounding box to the current shape.
#'
#' @description Take a shapefile as a sf object and clips boxes from a Halton frame around it. Size of those boxes is chosen
#' by choosing J, the number of base 2,3 powers to subdivide. Intended for internal use but can be useful in other
#' context.
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
#'
#' @param shp shape as spatial features object to wrap into Halton frame.
#' @param bb Master Sample bounding box.
#' @param base Co-prime base for BAS, do not change from 2,3.
#' @param J Definition for the number of grid cells of Halton frame.
#' @param projstring Projection that the master sample is in, can be passed as part of the bounding box.
#' @param rotate Boolean if you want to rotate shape before exporting.
#'
#' @return rotated sf spatial object.
#'
#' @examples
#' \dontrun{
#' bb <- uc511::getBB()
#' data(NS_bioregion)
#' haltonBoxes <- uc511::shape2Frame(shp = NS_bioregion, J = c(6,4), bb = bb)
#' }
#'
#' @export
shape2Frame <- function(shp, bb = NULL, base = c(2,3), J = c(2,2), projstring = NULL, rotate = FALSE){

  if( !is.null( bb)){
    bb.bounds <- sf::st_bbox(bb)
    scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
    shift.bas <- bb.bounds[1:2]
    theta <- base::attr(bb, "rotation")
    cntrd <- base::attr(bb, "centroid")
    projstring <- sf::st_crs(bb)$proj4string
  }else{ return("uc511(shape2Frame) Define Bounding Box Please.")}

  if(base::is.null(projstring)) {
    projstring <- uc511::getProj()
    base::message("uc511(shape2Frame) Assuming Projection\n")
  }
  if(sf::st_crs(shp) != sf::st_crs(projstring)) shp <- sf::st_transform(shp, projstring)

  #Stretch bounding box to Halton Frame Size:
  shp2 <- uc511::rotate.shp(shp, bb, back = FALSE)
  bb2 <- sf::st_bbox(shp)
  xy <- (bb2 - shift.bas[c(1,2,1,2)])/scale.bas[c(1,2,1,2)]
  lx <- base::floor(xy[1] / (1/base[1]^J[1]))/(base[1]^J[1])
  ly <- base::floor(xy[2] / (1/base[2]^J[2]))/(base[2]^J[2])
  ux <- base::ceiling(xy[3] /(1/base[1]^J[1]))/(base[1]^J[1])
  uy <- base::ceiling(xy[4] /(1/base[2]^J[2]))/(base[2]^J[2])
  nx <- (ux-lx)*base[1]^J[1]
  ny <- (uy-ly)*base[2]^J[2]

  bb.new <- c(lx,ly, ux, uy)*scale.bas[c(1,2,1,2)] + shift.bas[c(1,2,1,2)]
  halt.frame <- raster::raster(raster::extent(base::matrix( bb.new , 2, 2)), nrow=ny, ncol=nx)
  raster::projection(halt.frame) <- projstring
  halt.poly <- raster::rasterToPolygons(halt.frame)
  if(rotate) return(uc511::rotate.shp(sf::st_as_sf(halt.poly), bb))

  return(sf::st_as_sfc(halt.poly))
}
