###########################################################################
## sea.R
## Fits a matern covariance function so spatial data using Integrated Nested Laplace Approximation (INLA)
## Author: Colin Millar
## This code comes with no warranty or guarantee of any kind.
###########################################################################

setwd("n:/colin/2011/11_11_soap_whiting")

#source("http://www.math.ntnu.no/inla/givemeINLA.R")
#install.packages(c("fields","rgl"))

require(INLA)
require(rgl)
require(fields)
require(grid)
source("scripts/data-functions.r")

########################
## fit function inputs: 
########################

## data preparation - a data frame directly from DATRAS with extra columns for observations (e.g. log recruitment)

ibts.fnames <- paste("data", c("NSibts.csv", "WCibts.csv"), sep="/")
ibts <- do.call(rbind, lapply(ibts.fnames, read.table, sep=",", header = TRUE, stringsAsFactors = FALSE))

# get rectuits
ibts $ rec <- ibts $ Age_1

# work with logged data, substitute zeros with minumum observed...
ibts $ log.rec <- with(ibts, log(ifelse(rec == 0, min(rec[rec >0]), rec)))

## The boundary - points should be organised clockwise
#               - 1st polygon is outer boundary
#               - any other polygons are 'islands'
#               - care should be taken so that these do not overlap

load("data/limits.rdata")

# make 3D
limits3D <- lapply(limits, lonlat3D)

mesh <- create_surface(limits, max.edge = 1) 
mesh3D <- create_surface(limits3D, max.edge = 0.015) 
plot( mesh )
plot( mesh3D )

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


##############################
# get subset of data and plot
##############################

Data <- 
  with(subset(ibts, Year == 2000 & Quarter == 1 & grepl("Mer", ibts $ Species)), 
    list(loc = cbind(x = ShootLon, y = ShootLat), y = log.rec))

# make 3D
Data3D <- Data
Data3D $ loc <- lonlat3D(Data $ loc)	
	
plot.data(Data3D $ loc, mesh = mesh3D, land = land3D, depths = contours3D, grd = grd3D)
# with no mesh
plot.data(Data3D $ loc, land = land3D, depths = contours3D, grd = grd3D)
# with no mesh but with boundary
plot.data(Data3D $ loc, land = land3D, depths = contours3D, limits = limits3D, grd = grd3D)


##########
## fit
##########

fit <- fit_surface(Data3D, mesh3D)


##########
## predict
##########

# projection for prediction
proj3D <- inla.mesh.projector.cm(mesh3D, dims=c(315,216)*2, projection='longlat', xlim = c(-9.5, 8.5), ylim = c(50.5, 62))

# predict from model
z <- exp( inla.mesh.project(proj3D, fit $ mean))


##########
## plot
##########

# simple plot
image(z)

# set sensible break and contour levels
breaks <- c(seq(0, .4, length = 180), seq(.5, 1, by = .1)) * 50000 
clevels = c(10, 100, 1000, 2500, 5000, 10000, 20000)

# pretty plot
plot.fit(z, proj, colfun = cm.cols2, breaks = breaks, clevels = clevels, land = land, limits = limits)

# plot standard errors (on log scale)
plot.fit(inla.mesh.project(proj, fit $ sd), proj)

# plot .025 quantiles
plot.fit(exp( inla.mesh.project(proj, fit[["0.025quant"]])), proj, clevels = clevels)

## tada!

