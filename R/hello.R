# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

#BASMastersample()

#- SRS:                    SRS
#- BAS:                    BAS (equal, unequal, seed point, panels?)
#- HaltonFrames:           Halton frames (for discretizing a continuous resource)
#- HIP:                    HIP (equal probability)
#- SpacialBalanceMeasures: Spatial balance measures (Voronoi and Modified Moran’s I)
#- VarianceEstimators:     Variance estimators (local mean, nearest neighbour and quasi bootstrap)
#- DAS:                    DAS
#- SpatialStratification:  Spatial stratification algorithm?
#- BASMasterSample:        BASMastersample

# to create a vignette run the following (just once!):
# usethis::use_vignette("introduction")

#devtools::install_github("paul-vdb/DFO-master-sample")


#library(BASMasterSample)
#library(sf)
#library(sp)
#library(raster)

#data(Fed_MPAs_clipped)
#library(uc511)
#library(profvis)
#profvis(
#smp <- xmasterSample(Fed_MPAs_clipped, N = 1000)
#)

# boxes <- ifelse(boxes < boxInit, B + (boxes - boxInit), boxes - boxInit)
#boxes = c(2, 29, 34, 25, 22, 13,  1)
#boxInit = c(35)
#B = 36
#B + (boxes - boxInit)


#plot(smp)
#plot(Fed_MPAs_clipped, max.plot = 40)

#library(rgdal)
#xx <- readOGR(dsn=Fed_MPAs_clipped[1], layer = "Area1")
#library(raster)
#extent(Fed_MPAs_clipped)
#class(Fed_MPAs_clipped)
#crs(Fed_MPAs_clipped)
