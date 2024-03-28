# getHaltonPointsFromExpandableGrid.R

#' @export
getHaltonPointsFromExpandableGrid <- function(shapefile,
                                              N,
                                              J,
                                              bases,
                                              seeds,
                                              crs,
                                              verbose = FALSE){
  # initialise variables.
  i <- 0
  pts_in_intersection <- 0
  # run loop until we have sufficient points...
  while (pts_in_intersection <= N){

    # get the frame.
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

    if(verbose){
      msg <- "uc511(getHaltonPointsFromExpandableGrid) Points in intersection: %s."
      msgs <- sprintf(msg, pts_in_intersection)
      base::message(msgs)
    }

    # expand the Halton frame.
    i <- i + 1
  } # end while

  if(verbose){
    # display some statistics before returning results.
    msg <- "uc511(getHaltonPointsFromExpandableGrid) %s samples found in %s iterations, using J1=%s and J2=%s."
    msgs <- sprintf(msg, pts_in_intersection, i, J[1]+i-1, J[2]+i-1)
    base::message(msgs)
  }

  result <- base::list(i                   = i,
                       diff_               = diff_,
                       pts.shp             = pts.shp,
                       bb.new              = bb.new,
                       seeds               = seeds,
                       pts_in_intersection = pts_in_intersection)
  return(result)
}

