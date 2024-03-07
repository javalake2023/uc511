# SRS.R

#' @name SRS
#'
#' @title Simple random sampling.
#'
#' @description This function invokes base::sample() to draw a random sample using
#' a user specified random seed.
#'
#' @details This function was written by Phil Davies.
#'
#' @param seed The random seed to be used to draw the current sample.
#' @param total_rows The total number of rows in our input dataset.
#' @param sample_size The number of rows wanted in our random sample.
#'
#' @return A random sample.

#' @examples
#' \dontrun{
#' # Create a random sample with a seed of 99.
#' uc511::SRS(seed = 99, total_rows = 100, sample_size = 20)
#' # Create a random sample with a seed of 42.
#' uc511::SRS(seed = 42, total_rows = 100, sample_size = 20)
#' # Create a random sample with a seed of 99.
#' uc511::SRS(seed = 99, total_rows = 100, sample_size = 25)
#' }
#'
#' @export
SRS <- function(seed = 42, total_rows = 0, sample_size = 0) {

  # validate our parameters.
  uc511::validate_parameters("seed", c(seed))
  uc511::validate_parameters("total_rows", c(total_rows))
  uc511::validate_parameters("sample_size", c(sample_size))
  # ensure sample_size < total_rows
  if(sample_size >= total_rows){
    stop("uc511(SRS) Parameter sample_size must be less than total_rows.")
  }

  samp_pts <- NULL
  base::set.seed(seed)
  samp_pts <- base::sample(x      = total_rows,
                          size    = sample_size,
                          replace = FALSE,
                          prob    = NULL)
  return (samp_pts)
}
