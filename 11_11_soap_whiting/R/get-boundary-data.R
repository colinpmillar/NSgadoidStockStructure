###########################################################################
## Author: Colin Millar
## This code comes with no warranty or guarantee of any kind.
###########################################################################

setwd("~/work/11_11_soap_whiting")

source("R/tools.R")

## The boundary - points should be organised clockwise
#               - 1st polygon is outer boundary
#               - any other polygons are 'islands'
#               - care should be taken so that these do not overlap

load("data/limits.rdata")

# make 3D
limits3D <- lapply(limits, lonlat3D)

mesh <- create_surface(limits, max.edge = 1) 
mesh3D <- create_surface(limits3D, max.edge = 0.015) 

## land and depth contours for plotting only
load("data/land.rdata")
load("data/contours.rdata")

# get lat lon grid
grd <- 
  rbind(
    expand.grid(x = c(seq(-10, 9, by = .01), NA), y = seq(51, 61, by = 1)),
    expand.grid(y = c(seq(51, 61, by = .01),NA), x = seq(-10, 9, by = 1))[2:1] )


# make 3D
land3D <- lapply(land, lonlat3D)
contours3D <- lapply(contours, lonlat3D)
grd3D <- lonlat3D(grd)

# projection for prediction
proj3D <- inla.mesh.projector.cm(mesh3D, dims=c(315,216)*2, projection='longlat', xlim = c(-9.5, 8.5), ylim = c(50.5, 62))


save(land3D, contours3D, grd3D, proj3D, mesh3D, limits3D, mesh,  file="run/boundary-data.RData")

