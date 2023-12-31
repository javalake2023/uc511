% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cppWhere2Start}
\alias{cppWhere2Start}
\title{Internal function to find the ordering of the first box according to the random seed.}
\usage{
cppWhere2Start(J, seeds, bases, boxes)
}
\arguments{
\item{J}{Definition for the number of grid cells of Halton frame.}

\item{seeds}{Master Sample random seed.}

\item{bases}{Co-prime bases should really always be 2,3}

\item{boxes}{ordering of boxes that have been clipped to be reordered according to the master sample seed.}
}
\value{
vector of reordered Halton indices.
}
\description{
This is a function to find which Halton Box the initial BAS point from the Master Sample falls into and
thus use it to order the remaining boxes based on the initial. It also helps us tracks
the master sample index as we skip boxes that have no resource.
}
\details{
This function was first written in R by Paul van Dam-Bates for the
package BASMasterSample. Subsequently it was written in C/C++ by Phil Davies.
}
