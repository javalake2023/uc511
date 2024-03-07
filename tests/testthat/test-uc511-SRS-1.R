# test-uc511-SRS-1.R


testthat::test_that("1. Test SRS.", {
  samp1 <- SRS(total_rows = 10, sample_size = 3)
  expect_equal(samp1, c(1, 5, 10))
})

testthat::test_that("2. Test SRS.", {
  samp2 <- SRS(seed = 69, total_rows = 10, sample_size = 3)
  expect_equal(samp2, c(1, 2, 8))
})

testthat::test_that("3. Test SRS.", {
  samp3 <- SRS(seed = 136, total_rows = 100, sample_size = 4)
  expect_equal(samp3, c(60, 9, 85, 34))
})
