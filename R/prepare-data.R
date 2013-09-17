

setwd("n:/colin/2011/11_11_soap_whiting")

	   
get.land.masses <- 
function()
{
  load("data/gmt4.rda")

  xy <- gmt4
  names(xy) <- c("x","y")

  change <- c(0, diff(as.numeric(is.na(xy $ x))))
  stop <- c(which(change == 1) - 1, nrow(xy))
  start <- c(1, which(change == -1))

  xy.list.full <- lapply(1:length(start), function(i) xy[start[i]:stop[i], ])
  lens <- sapply(xy.list.full, nrow)

  c(xy.list.full[lens > 250], xy.list.full[c(270,274,284,332,318,327,328,331,335,336,339,375,398,399,438)])
}

get.contours <-
function()
{
  depth <- read.table("data/northsea-depths.xyz")
  names(depth) <- c("x","y","depth")

  dxs <- unique(depth $ x)
  dys <- rev(unique(depth $ y))
  depthz <- matrix(depth $ depth, length(dxs), length(dys))			   
  depthz <- depthz[,length(dys):1]

  clines <- contourLines(dxs, dys, depthz, levels = -200)
  clines <- lapply(clines, function(x) with(x, data.frame(x=x, y=y)))
  contours <- clines[c(1,47)]
  clines <- contourLines(dxs, dys, depthz, levels = -10)
  clines <- lapply(clines, function(x) with(x, data.frame(x=x, y=y)))[sapply(clines, function(x) length(x$x)) > 100]
  c(contours, clines[-c(4,5,12,14,15)]) 
}

 
# get coastline data
land <- get.land.masses()
contours <- get.contours()

save(land, file = "data/land.rdata")
save(contours, file = "data/contours.rdat")


# get boundaries
bndry = read.table("data/new-bndry.dat", header = TRUE)
bndry = as.matrix(rbind(bndry, bndry[1,]))
hebs = read.table("data/hebs.dat", header = TRUE)
hebs = as.matrix(rbind(hebs, hebs[1,]))
shet = read.table("data/shetland.dat", header = TRUE)
shet = as.matrix(rbind(shet, shet[1,]))
ork = read.table("data/orkney.dat", header = TRUE)
ork = as.matrix(rbind(ork, ork[1,]))

limits <- list(bndry = bndry, hebs = hebs, shet = shet, ork = ork)

save(limits, file = "data/limits.rdata")

# useful code for creating boundaries

if (0)
{
bndry <- data.frame(x = numeric(0), y = numeric(0))
for (i in 1:10000)
{
plot.data()
lines(bndry, lwd = 2)
bndry <- rbind(bndry, as.data.frame(locator()))
}
write.table(bndry, "data/new.bndry.dat", row.names = FALSE)
}