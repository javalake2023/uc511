# Validate uc511 functions, features and parameter validation.

testthat::test_that("1. Verify internal function functions correctly.", {
  sf_object <- sf::st_read(system.file("shape/nc.shp", package="sf"))
  expect_error(uc511::getOverlappedPoints(sf_object, 1),
               "uc511(getOverlappedPoints) Simple file object does not contain a feature named panel_id.", fixed=TRUE)
})
