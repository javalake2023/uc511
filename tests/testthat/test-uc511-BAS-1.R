# Validate BAS functions, features and parameter validation.

testthat::test_that("1. Verify panels= and panel_overlap= parm length checking.", {
  n_panels <- c(20, 20, 20, 100)
  n_panel_overlap <- c(0, 4, 5)
  shp.cant <- NULL
  bb <- NULL
  expect_error(uc511::BAS(shp.cant,
                          panels = n_panels,
                          panel_overlap = n_panel_overlap,
                          bb = bb), "uc511(ValidatePanelDesign) length of panels [4] must match length of panel_overlap [3].", fixed=TRUE)
})

testthat::test_that("2. Verify panels= and panel_overlap= parm length checking.", {
  n_panels <- c(20, 20, 20)
  n_panel_overlap <- c(0, 4, 5, 6)
  shp.cant <- NULL
  bb <- NULL
  expect_error(uc511::BAS(shp.cant,
                          panels = n_panels,
                          panel_overlap = n_panel_overlap,
                          bb = bb), "uc511(ValidatePanelDesign) length of panels [3] must match length of panel_overlap [4].", fixed=TRUE)
})



