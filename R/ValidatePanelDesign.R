#' @name ValidatePanelDesign
#'
#' @title Validate the panels and panel_overlap parameters.
#'
#' @description This function is used to validate the panels and panel_overlap
#' parameters. The panel_design flag is set TRUE when the panels and/or panel_overlap
#' parameters are not NULL. This is an internal only function.
#'
#' @details This function was written by Phil Davies.
#'
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
#' @param n The number of samples required. Used when panels and panel_overlap are NULL.
#'
#' @return Returns a list of the following variables: n
#'                                                    panel_design
#'                                                    number_panels
#'                                                    panel_overlap

ValidatePanelDesign <- function(panels, panel_overlap, n){
  #
  # default is not a panel design. Will be set to true if either of the
  # panels or panel_overlap parameters are not NULL.
  panel_design <- FALSE
  number_panels <- 0

  # verify panels parameter, must be a list of numerics (if not null).
  # if panels not NULL, then we will ignore the n parameter.
  if(!is.null(panels)){
    validate_parameters("panels", panels)
    n <- base::sum(panels)
    panel_design <- TRUE
    number_panels <- length(panels)
  }

  # verify panels_overlap parameter, must be a list of numerics (if not null).
  if(!is.null(panel_overlap)){
    validate_parameters("panel_overlap", panel_overlap)
    if(!panel_design){
      stop("uc511(ValidatePanelDesign) panels parameter must be specified when panel_overlap specified.")
    }
    if(length(panels) != length(panel_overlap)){
      msg <- "uc511(ValidatePanelDesign) length of panels [%s] must match length of panel_overlap [%s]."
      msgs <- sprintf(msg, length(panels), length(panel_overlap))
      stop(msgs)
    }
    panel_design <- TRUE
    # force zero for panel 1.
    panel_overlap[1] <- 0
  }
  result <- base::list(panel_design  = panel_design,
                       number_panels = number_panels,
                       panel_overlap = panel_overlap,
                       n             = n)
  return(result)
}

#' @name PanelDesignAssignPanelids
#'
#' @title Assign panel ids to the samples.
#'
#' @description This function assigns panel id's to each sample based on values in the
#' panels and panel_overlap parameters. This is an internal only function.
#'
#' @details This function was written by Phil Davies.
#'
#' @param smp A shapefile for the region under study.
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
#' @param panel_design A flag, when TRUE, indicates that we are performing a panel design and
#' the parameters used are specified in the panels and panel_overlap parameters.
#' @param number_panels The number of sample panels required.
#'
#' @return Returns a list of the following variables: sample, which has the appropriate panel id
#' assigned determined by the panels and panel_overlap parameters.
#'

PanelDesignAssignPanelids <- function(smp, panels, panel_overlap, panel_design, number_panels){
  # if(panel_design) then assign panel id's to smp.
  # need to distinguish if panel_overlap is required.
  if(is.null(panel_overlap) & panel_design){
    tmp <- NULL
    # assign panel id's to sample points. smp$panel_id.
    for(i in 1:number_panels){
      tmp <- c(tmp, rep(i, panels[i]))
    }
    smp$panel_id <- tmp
  }
  # if(panel_overlap) is not null and panel_design is TRUE
  if(!is.null(panel_overlap) & panel_design){
    message("panel_overlap")
    # need to create the panel_id column
    # Initialize variables
    smp$panel_id <- 0
    panelid <- 1
    start_index <- 1

    for(i in 1:length(panels)){
      start_index <- start_index - panel_overlap[i]
      if(i == 1){
        for(j in 1:panels[i]){
          smp$panel_id[start_index] <- panelid
          start_index <- start_index + 1
        }
      } else {
        for(j in 1:panels[i]){
          #print(unlist(df$panelid[start_index]))
          #if(df$panelid[start_index] == 0){
          #  print("just a zero.")
          #}
          smp$panel_id[start_index] <- list(c(smp$panel_id[start_index], panelid))
          start_index <- start_index + 1
        }
      }
      panelid <- panelid + 1
    }
  }
  result <- base::list(sample = smp)
  return(result)
}

