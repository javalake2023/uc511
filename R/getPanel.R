# getPanel.R

#' @name contains_feature
#'
#' @title Check if the sf object contains a specified feature.
#'
#' @description Used to check if a simple file object contains a feature. This
#' is an internal only function.
#'
#' @details This function was written by Phil Davies.
#'
#' @param sf_object Simple file object that we want to verify if it contains
#' a feature called feature_name.
#' @param feature_name The feature name we want to find in the simple file
#' object sf_object.
#'
#' @return Returns TRUE if the simple file object sf_object contains the feature
#' feature_name. Otherwise FALSE is returned.
#'
#' @examples
#' \dontrun{
#' # Check if a feature is available in a shapefile.
#' containsFeature <- contains_feature(shapefile, "NAME")
#' }
#'
contains_feature <- function(sf_object, feature_name) {
  # Drop the geometry column
  df <- sf::st_drop_geometry(sf_object)

  # Check if the feature exists in the data frame
  feature_exists <- feature_name %in% base::names(df)

  return(feature_exists)
}


#' @name getPanel
#'
#' @title Select points from a polygon using a BAS Master Sample.
#'
#' @description This is the main function for selecting sites using the BAS master
#' sample. It assumes that you have already defined the master sample using the
#' BoundingBox() function or will be selecting a marine master sample site in BC.
#'
#' @details This function was written by Phil Davies.
#'
#' @param shapefile Shape file as a polygon (sp or sf) containing a sample that
#' contains a feature column named panel_id.
#' @param panelid The overlapped panel in the shapefile shp the user wants
#' sample points from.
#'
#' @return The sample for the specified panel.

#' @examples
#' \dontrun{
#' # Get all the sample from panel 1.
#' panelid <- 1
#' panel_1 <- uc511::getPanel(shp, panelid)
#' }
#'
#' @export
getPanel <- function(shapefile, panelid){
  # validate our parameters. Ensure panelid is numeric and has a value greater than zero.
  uc511::validate_parameters("panelid", c(panelid))

  # Usage: Check if the sf object contains a feature named "panel_id"
  feature_name <- "panel_id"
  exists <- contains_feature(shapefile, feature_name)
  if(!exists){
    base::stop("uc511(getPanel) Simple file object does not contain a feature named panel_id.")
  }

  indx <- base::c()
  for(k in 1:base::length(shapefile$panel_id)){
    tmp <- shapefile$panel_id[k]
    if(base::any(panelid %in% base::unlist(tmp))){
      indx <- base::c(indx, k)
    }
  }
  smp <- shapefile[indx,]
  result <- base::list(sample  = smp,
                       indices = indx)
  return(result)
}
