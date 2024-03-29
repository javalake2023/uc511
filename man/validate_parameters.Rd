% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ValidateParameters.R
\name{validate_parameters}
\alias{validate_parameters}
\title{Validate uc511 function parameters.}
\usage{
validate_parameters(parm, parm_value)
}
\arguments{
\item{parm}{The parameter to be validated.}

\item{parm_value}{The value of the parameter to be validated. Must be defined as a list.}
}
\value{
Always returns TRUE indicating that the parameter was parsed successfully. If
a parameter fails validation further execution is terminated using the STOP function.
}
\description{
This function is used to validate parameters passed to
all available uc511 functions.
}
\details{
This function was written by Phil Davies.
}
