# Validate HIP functions.

test_that("Compare cppRSHalton() and HaltonPts().", {
  chp <- uc511::cppRSHalton(n = 100000)
  hp <- HaltonPts(n = 100000)
  expect_equal(chp[,2:3], hp)
})
