###########################################################################
## Fits a matern covariance function so spatial data using Integrated Nested Laplace Approximation (INLA)
## Author: Colin Millar
## This code comes with no warranty or guarantee of any kind.
###########################################################################

setwd("~/work/11_11_soap_whiting")

#source("http://www.math.ntnu.no/inla/givemeINLA.R")
#install.packages(c("fields","rgl"))

require(INLA)
require(fields)
source("R/tools.R")

########################
## fit function inputs: 
########################

load("run/boundary-data.RData")
load("run/whiting-ibts.RData")

##############################
# fit and predict function
##############################

doone <- function (subyear, what) {

  # data required for a fit
  Data <- 
    with(subset(ibts, year == subyear), 
         list(loc = cbind(x = shootlon, y = shootlat), rec = rec, ssb = ssb))

  # make 3D and get response
  Data3D <- list(loc = lonlat3D(Data $ loc), y = Data[[what]])

  # work with logged data, substitute zeros with minumum observed... can always do better here...
  Data3D $ y <- with(Data3D, log(ifelse(y == 0, min(y[y >0]), y)))
	
  ## fit
  cat("fitting", what, subyear, "\n")
  fit <- fit_surface(Data3D, mesh3D)

  ## predict
  exp( inla.mesh.project(proj3D, fit $ mean))
}

years <- c(0, 5, 10) + 1990

rec.fits <- lapply(years + 1, doone, what = "rec")
ssb.fits <- lapply(years, doone, what = "ssb")

save(rec.fits, ssb.fits, years, file = "run/output.RData")


