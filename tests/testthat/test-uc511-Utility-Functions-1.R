# Validate uc511 functions, features and parameter validation.

testthat::test_that("1. Verify internal function functions correctly.", {
  sf_object <- sf::st_read(system.file("shape/nc.shp", package="sf"))
  expect_error(uc511::getOverlappedPoints(sf_object, 1),
               "uc511(getOverlappedPoints) Simple file object does not contain a feature named panel_id.", fixed=TRUE)
})

# Validate and compare Halton point generators.

testthat::test_that("2. Compare cppRSHalton() and cppRSHalton_br() same n, seeds.", {
  chp <- uc511::cppRSHalton(n = 1000, seeds = c(123, 456))
  chpbr <- uc511::cppRSHalton_br(n = 1000, seeds = c(123, 456))
  expect_equal(chp[,2:3], chpbr$pts)
})

testthat::test_that("3. Compare cppBASpts() and cppRSHalton_br() same n, bases.", {
  chp <- uc511::cppBASpts(n = 1000, bases = c(2, 3))
  chpbr <- uc511::cppRSHalton_br(n = 1000, bases = c(2, 3), seeds = chp$seeds)
  expect_equal(chp$pts, chpbr$pts)
})

