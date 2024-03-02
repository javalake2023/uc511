# Validate uc511 functions, features and parameter validation.

testthat::test_that("1. Verify validate_parameters parses parmeters correctly.", {
  expect_error(uc511::validate_parameters("testparm", "A"),
               "uc511(validate_parameters) Parameter testparm must contain all numeric values.", fixed=TRUE)
})

testthat::test_that("2. Verify validate_parameters parses parmeters correctly.", {
  expect_error(uc511::validate_parameters("J", 99),
               "uc511(validate_parameters) Parameter J must be a list of length 2.", fixed=TRUE)
})

testthat::test_that("3. Verify validate_parameters parses parmeters correctly.", {
  expect_error(uc511::validate_parameters("hipIterations", 1),
               "uc511(validate_parameters) Parameter hipIterations values less than two are not supported.", fixed=TRUE)
})

testthat::test_that("4. Verify validate_parameters parses parmeters correctly.", {
  expect_error(uc511::validate_parameters("hipIterations", 14),
               "uc511(validate_parameters) Parameter hipIterations values greater than 13 are not supported.", fixed=TRUE)
})

testthat::test_that("5. Verify validate_parameters parses parmeters correctly.", {
  expect_error(uc511::validate_parameters("panelid", 0),
               "uc511(validate_parameters) Parameter panelid must have a value greater than zero.", fixed=TRUE)
})



