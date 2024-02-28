#' @name getOverlappedPoints
#'
#' @title Select points from a polygon using a BAS Master Sample.
#'
#' @description This is the main function for selecting sites using the BAS master
#' sample. It assumes that you have already defined the master sample using the
#' BoundingBox() function or will be selecting a marine master sample site in BC.
#'
#' @details This function was written by Phil Davies.
#'
#' @param shp Shape file as a polygon (sp or sf) containing a BAS sample that
#' contains panelid's.
#' @param panelid The overlapped panel in the shapefile shp the user wants
#' sample points from.
#'
#' @return The sample for the specified panel.

#' @examples
#' \dontrun{
#' #
#' }
#'
#' @export
getOverlappedPoints <- function(shapefile, panelid){
  # validate our parameters. Ensure panelid is numeric and has a value greater than zero.
  uc511::validate_parameters("panelid", c(panelid))

  indx <- c()
  for(k in 1:length(shapefile$panel_id)){
    tmp <- shapefile$panel_id[k]
    if(any(panelid %in% unlist(tmp))){
      indx <- c(indx, k)
    }
  }
  smp <- shapefile[indx,]
  result <- base::list(sample  = smp,
                       indices = indx)
  return(result)
}
