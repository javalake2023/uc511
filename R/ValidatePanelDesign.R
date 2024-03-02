#


ValidatePanelDesign <- function(panels, panel_overlap, n){
  #
  # default is not a panel design. Will be set to true if either of the
  # panels or panel_overlap parameters are not NULL.
  panel_design <- FALSE

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


PanelDesignAssignPanelids <- function(smp, panels, panel_overlap, panel_design, number_panels){
  # if(panel_design) then assign panel id's to smp.
  # need to distinguish if panel_overlap is required.
  if(panel_design & is.null(panel_overlap)){
    tmp <- NULL
    # assign panel id's to sample points. smp$panel_id.
    for(i in 1:number_panels){
      tmp <- c(tmp, rep(i, panels[i]))
    }
    smp$panel_id <- tmp
  }
  # if(panel_overlap) is not null and panel_design is TRUE
  if(!is.null(panel_overlap) & panel_design){
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

