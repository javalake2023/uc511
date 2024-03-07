# Validate HIP functions.

testthat::test_that("1. Place holder until we have HIP stuff to test.", {
  chp <- uc511::cppRSHalton(n = 1000, seeds = c(123, 456))
  chpbr <- uc511::cppRSHalton_br(n = 1000, seeds = c(123, 456))
  expect_equal(chp[,2:3], chpbr$pts)
})

