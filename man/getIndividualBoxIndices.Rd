% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MasterSample.R
\name{getIndividualBoxIndices}
\alias{getIndividualBoxIndices}
\title{Get the Halton Box Index for a point.}
\usage{
getIndividualBoxIndices(pts, J = NULL, bb, size = 100)
}
\arguments{
\item{pts}{Shape file as a polygon (sp or sf) to generate the boxes for.}

\item{J}{Definition for the number of grid cells of Halton frame.}

\item{bb}{Bounding box which defines the Master Sample. Default is BC Marine Master Sample.}

\item{size}{Physical target size of Halton boxes (square ish) if J is NULL.}
}
\value{
Return the Halton Index for all "pts" that are passed to this function.
}
\description{
If you pass discrete points to this function it will return the Halton Index according to the
Master Sample as defined by the bounding box. If the points are closer together than the defined box
size then they are treated as the same sample unit. We do not expect this to be spatially balanced and do
not recommend using this function alone to select a sample.
}
\details{
This function was first written by Paul van Dam-Bates for the
package BASMasterSample.
}
\examples{
\dontrun{
library(bcmaps)
# If you haven't already installed this:
install.packages('bcmapsdata', repos='https://bcgov.github.io/drat/')
library(bcmapsdata)
cities <- bcmaps::get_layer("bc_cities")
bb <- uc511::buildMS(cities, d = 2, FALSE)
# For visibility will make boxes 10 km
cities.ord <- uc511::getIndividualBoxIndices(pts = cities, bb = bb, size = 100)
plot(st_geometry(cities), pch = 20)
plot(st_geometry(cities.ord[rank(cities.ord$HaltonIndex) < 15,]), add= TRUE, col = "red", cex = 1.5)
}

}
