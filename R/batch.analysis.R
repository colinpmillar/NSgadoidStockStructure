
setwd("~/work/11_11_soap_whiting")

#source("http://www.math.ntnu.no/inla/givemeINLA.R")
#install.packages(c("fields","rgl"))

require(INLA)
require(rgl)
require(fields)
require(grid)
source("scripts/data-functions.r")


## out directory
outdir <- "figures"

## data preparation - a data frame directly from DATRAS with extra columns for observations (e.g. log recruitment)

ibts.fnames <- paste("data", c("NSibts.csv", "WCibts.csv"), sep="/")
ibts <- do.call(rbind, lapply(ibts.fnames, read.table, sep=",", header = TRUE, stringsAsFactors = FALSE))

# get rectuits and SSB
stats <- list(rec = c(1, 0,     0,     0,     0,     0),
              #ssb = c(0, 0.177, 0.214, 0.248, 0.287, 0.299)) # WC medians 
			  ssb = c(0, 0.179, 0.239, 0.285, 0.310, 0.348)) # NS medians
			  
## apply appropriate mats and weights
weight <- expand.grid(Species = unique(ibts $ Species), 
                      Survey = unique(ibts $ Survey), stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
mat <- expand.grid(Species = unique(ibts $ Species), 
                      Survey = unique(ibts $ Survey), stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

weight[paste("Age_", 0:10, sep="")] <- rbind(c(0   ,0    ,0.179,0.239,0.285,0.310,0.348,NA  ,NA  ,NA  ,NA   ), # NS Whi
                                             c(0   ,0.36 ,0.82 ,2.01 ,3.79 ,5.90 ,7.90 ,NA  ,NA  ,NA  ,NA   ), # NS Cod
										     c(0   ,0    ,0    ,0    ,0    ,0    ,0    ,0   ,0   ,0   ,0    ), # NS Had
					   					     c(0   ,0    ,0.177,0.214,0.248,0.287,0.299,NA  ,NA  ,NA  ,NA   ), # WC Whi
										     c(0   ,0.35 ,1.12 ,2.55 ,4.47 ,6.29 ,7.87 ,NA  ,NA  ,NA  ,NA   ), # WC Cod
										     c(0   ,0    ,0    ,0    ,0    ,0    ,0    ,0   ,0   ,0   ,0    )) # WC Had

mat[paste("Age_", 0:10, sep="")] <- rbind(c(0  ,0   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ), # NS Whi
                                          c(0  ,0.01,0.05,0.23,0.62,0.86,1   ,1   ,1   ,1   ,1   ), # NS Cod
										  c(0  ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ), # NS Had
										  c(0  ,0   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ), # WC Whi
										  c(0  ,0   ,0.52,0.86,1   ,1   ,1   ,1   ,1   ,1   ,1   ), # WC Cod
										  c(0  ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   )) # WC Had
										  
ibts $ ID <- 1:nrow(ibts)
weight <- merge(ibts[c("Species","Survey","ID")], weight, all = TRUE)
weight <- weight[order(weight $ ID),] 
mat <- merge(ibts[c("Species","Survey","ID")], mat, all = TRUE)
mat <- mat[order(mat $ ID),] 
					
cnames <- paste("Age_",0:10,sep="")
					
ibts $ ssb <- rowSums(weight[cnames] * mat[cnames] * ibts[cnames], na.rm = TRUE)
							

# my personal preference
names(ibts) <- tolower(names(ibts))
ibts <- ibts[c("year", "quarter", "species", "shootlat", "shootlon", "rec", "ssb")]


## The boundary - points should be organised clockwise
#               - 1st polygon is outer boundary
#               - any other polygons are 'islands'
#               - care should be taken so that these do not overlap

# get boundaries
bndry = read.table("data/westcoast-boundary.dat", header = TRUE)
bndry = as.matrix(rbind(bndry, bndry[1,]))
hebs = read.table("data/hebs.dat", header = TRUE)
hebs = as.matrix(rbind(hebs, hebs[1,]))
limits <- list(bndry = bndry, hebs = hebs)

mesh <- create_surface(limits, max.edge = .5) 

plot( mesh )

# projection for prediction
proj <- inla.mesh.projector(mesh, dims = c(315, 216) * 2)

## land and depth contours for plotting only
load("data/land.rdata")
load("data/contours.rdata")

ibts <- subset(ibts, shootlon < -4)

ibtsloc <- as.matrix(ibts[ibts $ quarter %in% c(1,4),c("shootlon","shootlat")])

# with no mesh but with boundary
plot.data(ibtsloc, land = land, depths = contours, limits = limits, mesh = mesh)



#############
## do batch
############

sp <- "morhua"
qtr <- 4
yr <- 1997
what <- "ssb"

out <- list()

for (sp in c("morhua")) {
for (what in c("ssb","rec")){
for (qtr in c(1, 4)) {
for (yr in 1986:2011) {

if (qtr == 4 & yr %in% c(1986:1996, 2010:2011)) next

if (qtr == 1 & yr %in% c(1988,1994) & what == "rec" & sp == "morhua") next

# data required for a fit
Data <- 
  with(subset(ibts, year == yr & quarter == qtr & grepl(sp, ibts $ species)), 
    list(loc = cbind(x = shootlon, y = shootlat), rec = rec, ssb = ssb))
Data $ y <- Data[[what]]

# work with logged data, substitute zeros with minumum observed...
Data $ y <- with(Data, log(ifelse(y == 0, min(y[y >0])/2, y)))
	
##########
## fit
##########

cat("qtr ==",qtr,"& yr ==",yr,"& what ==",  what, "& sp ==", sp, "\n"); flush.console()
fit <- fit_surface(Data, mesh)

##########
## predict
##########

# predict from model
z <- exp( inla.mesh.project(proj, fit $ mean))


##########
## plot
##########

# simple plot
fname <- paste("fit.",sp,".",what,".",yr,".Q",qtr,sep="")

plot.fit(z, proj, land = land, limits = limits)
points(Data $ loc, pch = 16, cex = .2)

#plot.fit(z, proj, land = land, limits = limits, fname = paste(outdir,fname,sep=""))

out[[fname]] <- z

}}}}

zsq1 <- out[grepl("ssb", names(out)) & grepl("Q1", names(out))]
zsq4 <- out[grepl("ssb", names(out)) & grepl("Q4", names(out))]

zrq1 <- out[grepl("rec", names(out)) & grepl("Q1", names(out))]
zrq4 <- out[grepl("rec", names(out)) & grepl("Q4", names(out))]

# set sensible break and contour levels
breaks <- c(seq(0, .4, length = 180), seq(.5, 1, by = .1)) * 800
#clevels = c(10, 100, 1000, 2500, 5000, 10000, 20000)

# pretty plot
plot.fit(zrq4[[1]], proj, colfun = cm.cols2, breaks = breaks, land = land, limits = limits)

# plot standard errors (on log scale)
#plot.fit(inla.mesh.project(proj, fit $ sd), proj)

# plot .025 quantiles
#plot.fit(exp( inla.mesh.project(proj, fit[["0.025quant"]])), proj, clevels = clevels)

## tada!



# recruit plots
breaks <- c(seq(0, .4, length = 180), seq(.5, 1, by = 5000)) 
cols <- cm.cols(length(breaks)-1)
clevels = c(10, 100, 1000, 2500, 5000, 10000, 20000) / 50000

if (0){
for (i in 1:length(fits))
{
  z <- inla.mesh.project(proj, fits[[i]] $ mean)
  tiff(paste("figures/recruits.", years[i], ".tiff",sep=""), 8, 8, res = 500, units = "in")
  plot.inla.fit(z, clevels, breaks, cols, FUN = exp)
  dev.off()  
}

for (i in 1:length(fits))
{
  z <- inla.mesh.project(proj, fits[[i]][["0.025quant"]])
  tiff(paste("figures/recruits.cil.", years[i], ".tiff",sep=""), 8, 8, res = 500, units = "in")
  plot.inla.fit(z, clevels, breaks, cols, FUN = exp)
  dev.off()  
}

for (i in 1:length(fits))
{
  z <- inla.mesh.project(proj, fits[[i]][["0.975quant"]])
  tiff(paste("figures/recruits.ciu.", years[i], ".tiff",sep=""), 8, 8, res = 500, units = "in")
  plot.inla.fit(z, clevels, breaks, cols, FUN = exp)
  dev.off()  
}

# sd (on log scale) of recruit plot
breaks <- seq(.5, 4.2, length = 200) 
cols <- cm.cols(length(breaks)-1)
clevels = pretty(breaks)

for (i in 1:length(fits))
{
  z <- inla.mesh.project(proj, fits[[i]] $ sd)
  tiff(paste("figures/recruits.sd.", years[i], ".tiff",sep=""), 8, 8, res = 500, units = "in")
  plot.inla.fit(z, clevels, breaks, cols)
  dev.off()  
}
}

# ssb plots
breaks <- c(seq(0, 200, length = 180), seq(300, 2000, by = 100)) 
cols <- cm.cols(length(breaks)-1)
clevels = c(10, 100, 1000, 2500, 5000, 10000, 20000)

if (0)
{
for (i in 1:length(fits))
{
  z <- inla.mesh.project(proj, fits[[i]] $ mean)
  tiff(paste("figures/ssb.", years[i], ".tiff",sep=""), 8, 8, res = 500, units = "in")
  plot.inla.fit(z, clevels, breaks, cols, FUN = exp)
  dev.off()  
}

for (i in 1:length(fits))
{
  z <- inla.mesh.project(proj, fits[[i]][["0.025quant"]])
  tiff(paste("figures/ssb.cil.", years[i], ".tiff",sep=""), 8, 8, res = 500, units = "in")
  plot.inla.fit(z, clevels, breaks, cols, FUN = exp)
  dev.off()  
}

for (i in 1:length(fits))
{
  z <- inla.mesh.project(proj, fits[[i]][["0.975quant"]])
  tiff(paste("figures/ssb.ciu.", years[i], ".tiff",sep=""), 8, 8, res = 500, units = "in")
  plot.inla.fit(z, clevels, breaks, cols, FUN = exp)
  dev.off()  
}

# sd (on log scale) of ssb
breaks <- seq(.5, 4.2, length = 200) 
cols <- cm.cols(length(breaks)-1)
clevels = pretty(breaks)

for (i in 1:length(fits))
{
  z <- inla.mesh.project(proj, fits[[i]] $ sd)
  tiff(paste("figures/ssb.sd.", years[i], ".tiff",sep=""), 8, 8, res = 500, units = "in")
  plot.inla.fit(z, clevels, breaks, cols)
  dev.off()  
}
}


