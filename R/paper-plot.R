###########################################################################
## Fits a matern covariance function so spatial data using Integrated Nested Laplace Approximation (INLA)
## Author: Colin Millar
## This code comes with no warranty or guarantee of any kind.
###########################################################################

setwd("~/work/11_11_soap_whiting")

#source("http://www.math.ntnu.no/inla/givemeINLA.R")
#install.packages(c("fields","rgl"))

require(INLA)
require(grid)
source("R/tools.R")

########################
## load data and fits
########################

load("run/boundary-data.RData")
load("run/whiting-ibts.RData")
load("run/output.RData")
load("data/land.rdata")
load("data/limits.rdata")

##########
## arrange data for plotting
##########

p.dat <- expand.grid(x = proj3D $ x, y = proj3D $ y, year = years)
p.dat $ rec <- unlist(rec.fits)
frec <- function(rec) rec / max(p.dat $ rec, na.rm = TRUE)
ifrec <- function(at) at * max(p.dat $ rec, na.rm = TRUE)
p.dat $ s.rec <- frec(p.dat $ rec)
p.dat $ ssb <- unlist(ssb.fits)
fssb <- function(ssb) ssb / max(p.dat $ ssb, na.rm = TRUE)
ifssb <- function(at) at * max(p.dat $ ssb, na.rm = TRUE)
p.dat $ s.ssb <- fssb(p.dat $ ssb)


scale.dat <- subset(p.dat, x > (mean(p.dat $ x) - 1) & x < (mean(p.dat $ x) + 1) & year == years[1])[c("x","y")]
fz <- function(y) (y - min(p.dat $ y)) / (max(p.dat $ y) - min(p.dat $ y)) * 1.1 - 0.05
ifz <- function(z) (z + 0.05) / 1.1 * (max(p.dat $ y) - min(p.dat $ y)) + min(p.dat $ y)
scale.dat $ z <- fz(scale.dat $ y)


ICES.areas <- list(x = c( -4,-4,-3,-3,-2, -2, 3, 3, 4, 4, 5, 5, 6, 6, 8, 8, 5, 5, 4,   4,   2, 2, 1,   1,  -1,-1,-2,-2,-3,-3,-2,  -2,  -4,  -4),
                   y = c(58.5,60,60,61,61,61.5,61.5,61,61,58.5,58.5,58,58,57.5,57.5,57,57,56,56,55.5,55.5,55,55,54.5,54.5,55,55,56,56,57,57,57.5,57.5,58.5))

split.line <- list(x = c(-4, -4, NA, 0.7389558, 7.8756676), y = c(58.60393, 59.9196, NA, 53.28031, 56.54514))

##########
## colour scheme
##########

cm.cols <- colorRampPalette(c("grey90","purple","cyan","yellow","red"))
cm.cols2 <- colorRampPalette(c("red", "pink","white"))
breaks <- c(seq(0, .4, length = 201), seq(.4, 1, length = 201)[-1]) 
cols <- c(cm.cols(200), cm.cols2(200))


##########
## plot elements
##########

between <- list(x = 0.2, y = 0.2)

legend.p <-
levelplot(s.rec ~ x * y | factor(year), data = subset(p.dat, year %in% years[1:2]),
          colorkey = FALSE, as.table = TRUE, 
		  xlab = "", ylab = "", layout = c(3,1), strip = FALSE,
		  between = between, 
		  par.settings = list(axis.line = list(col = "transparent")),
		  panel = panel.legend,
		  axis = axis.legend)
  
ssb.p <-
levelplot(s.ssb ~ x * y | factor(year), data = p.dat,
          colorkey = FALSE, as.table = TRUE, 
		  xlab = "", ylab = "", layout = c(3,1), strip = FALSE,
		  between = between,
		  par.settings = list(axis.line = list(col = "transparent")),
		  panel = panel.fun,  years = years,
		  axis = function(side, ...) axis.fun(side, plot.label = c("", "SSB"), ...))
		          
rec.p <-
levelplot(s.rec ~ x * y | factor(year), data = p.dat,
          colorkey = FALSE, as.table = TRUE, 
		  xlab = "", ylab = "", layout = c(3,1), strip = FALSE,
		  between = between,
		  par.settings = list(axis.line = list(col = "transparent")),
		  panel = panel.fun, years = years + 1,
		  axis = function(side, ...) axis.fun(side, plot.label = c("", "Recruitment (age 1)"), ...))



##########
## save plot
##########

scale <- 1.8
iheight <- 7.5 # inches
iwidth <- 7.5 # inches
dpi <- 300

png(file="fig/paperplot.png", width = iwidth*scale, height = iheight*scale, res = round(dpi/scale), units = "in")
  
# the figure
olap <- .06

plot(legend.p, position = c(0.2, .66 - olap, 1, 1), more = TRUE)
plot(ssb.p, position = c(0.2, .33 - olap/2, 1, .66 + olap/2), more = TRUE)
years <- years + 1
plot(rec.p, position = c(0.2, 0, 1, .33 + olap), more = FALSE)	  
years <- years - 1
dev.off()

system("convert fig/paperplot.png fig/paperplot.tiff")




