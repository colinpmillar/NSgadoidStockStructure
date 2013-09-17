
## plots!


## paper plot
pyears <- which(years %in% c(1986, 1991, 1996, 1998, 1999, 2000))

rec.z <- lapply(rec.fits[pyears], function(x) inla.mesh.project(proj, x $ mean))
ssb.z <- lapply(ssb.fits[pyears], function(x) inla.mesh.project(proj, x $ mean))

p.dat <- expand.grid(x = proj $ x, y = proj $ y, year = years[pyears])
p.dat $ rec <- exp(unlist(rec.z))
frec <- function(rec) rec / max(p.dat $ rec, na.rm = TRUE)
ifrec <- function(at) at * max(p.dat $ rec, na.rm = TRUE)
p.dat $ s.rec <- frec(p.dat $ rec)
p.dat $ ssb <- exp(unlist(ssb.z))
fssb <- function(ssb) ssb / max(p.dat $ ssb, na.rm = TRUE)
ifssb <- function(at) at * max(p.dat $ ssb, na.rm = TRUE)
p.dat $ s.ssb <- fssb(p.dat $ ssb)


scale.dat <- subset(p.dat, x > (mean(p.dat $ x) - 1) & x < (mean(p.dat $ x) + 1) & year == 1986)[c("x","y")]
fz <- function(y) (y - min(p.dat $ y)) / (max(p.dat $ y) - min(p.dat $ y)) * 1.1 - 0.05
ifz <- function(z) (z + 0.05) / 1.1 * (max(p.dat $ y) - min(p.dat $ y)) + min(p.dat $ y)
scale.dat $ z <- fz(scale.dat $ y)

ICES.areas <- list(x = c( -4,-4,-3,-3,-2, -2, 3, 3, 4, 4, 5, 5, 6, 6, 8, 8, 5, 5, 4,   4,   2, 2, 1,   1,  -1,-1,-2,-2,-3,-3,-2,  -2,  -4,  -4),
                   y = c(58.5,60,60,61,61,61.5,61.5,61,61,58.5,58.5,58,58,57.5,57.5,57,57,56,56,55.5,55.5,55,55,54.5,54.5,55,55,56,56,57,57,57.5,57.5,58.5))

split.line <- list(x = c(-4, -4, NA, 0.7389558, 7.8756676), y = c(58.60393, 59.9196, NA, 53.28031, 56.54514))

axis.fun <- function(side, plot.label, ...)
          {
            cprefix <- lattice:::lattice.getStatus("current.prefix")
            row <- lattice:::lattice.getStatus("current.focus.row", prefix = cprefix)
            column <- lattice:::lattice.getStatus("current.focus.column", prefix = cprefix)
		    switch(side,
              top = 
              {
			    if(row == 1 & column == 1)
                {				
			      grid.text(plot.label[1], x = unit(-.9, "npc"), y = unit(1.1, "npc"), just = "left") 
                  grid.text(plot.label[2], x = unit(-.1, "npc"), y = unit(1.1, "npc"), just = "right") 
				}
			  })
          }
		  
panel.fun <- function(x, y, z, subscripts, ...)
		  {
			cprefix <- lattice:::lattice.getStatus("current.prefix")
            row <- lattice:::lattice.getStatus("current.focus.row", prefix = cprefix)
            column <- lattice:::lattice.getStatus("current.focus.column", prefix = cprefix)
			tmp <- lapply(land, panel.polygon, col = "darkseagreen", border = FALSE)
			panel.levelplot(x, y, z, at = breaks, col.regions = cols, subscripts = subscripts)
			#panel.levelplot(x, y, z, at = clevels, contour = TRUE, region = FALSE, 
			#                subscripts = subscripts, col = grey(0.4))
			tmp <- lapply(segm_ns, function(seg) panel.polygon(seg $ loc, border = grey(0.3)))
			llines(split.line $ x, split.line $ y, col = "blue")
			grid.text(years[pyears[column + (row-1)*3]], x = unit(0.025, "npc"), y = unit(.95, "npc"), just = "left")
		  }

legend.p <-
levelplot(s.rec ~ x * y | factor(year), data = subset(p.dat, year < 1997),
          colorkey = FALSE, as.table = TRUE, 
		  xlab = "", ylab = "", layout = c(3,2), strip = FALSE,
		  between = list(x = 0.1, y = 0.2), 
		  par.settings = list(axis.line = list(col = "transparent")),
		  panel = function(x, y, z, subscripts, ...)
		  {
			cprefix <- lattice:::lattice.getStatus("current.prefix")
            column <- lattice:::lattice.getStatus("current.focus.column", prefix = cprefix)
			if (column == 2)
			{
			  tmp <- lapply(land, panel.polygon, col = "darkseagreen", border = FALSE)
			  t.sub = 1:nrow(mesh $ graph $ tv)
			  idx = cbind(mesh $ graph $ tv[t.sub, c(1:3, 1), drop = FALSE], NA)
              x = mesh$loc[t(idx), 1]
              y = mesh$loc[t(idx), 2]
              llines(x, y, col = grey(0.5), lwd = 1)
			  tmp <- lapply(segm_ns, function(seg) panel.polygon(seg $ loc, border = grey(0.3)))
			} else 
			if (column == 3)
			{
			  lpolygon(ICES.areas $ x, ICES.areas $ y)
			  tmp <- lapply(land, lpolygon, col = "darkseagreen", border = FALSE)
			  grid.text("IV", x = 2.72, y = 56.62, default.units = "native", just = "centre")
			} else 
			if (column == 1)
			{
			  # do scale...
			  panel.levelplot(scale.dat $ x, scale.dat $ y, scale.dat $ z, 
			                  at = breaks, col.regions = cols, subscripts = subscripts)
			  at <- seq(0, 1, by = 0.2)
			  rat <- pretty(ifrec(at))
			  sat <- pretty(ifssb(at))
			  grid.text(sprintf(rat, fmt = "%5.0f"), 
			             x = rep(min(scale.dat $ x) - .3, length(at)), 
			             y = ifz(frec(rat)), just = "right", default.units = "native", gp = gpar(cex = 0.7))
			  grid.text(sprintf(sat, fmt = "%4.0f"), 
			             x = rep(max(scale.dat $ x) + .3, length(at)), 
			             y = ifz(fssb(sat)), hjust = 0, default.units = "native", gp = gpar(cex = 0.7))
			}
		  },
		  axis = function(side, ...)
          {
            cprefix <- lattice:::lattice.getStatus("current.prefix")
            row <- lattice:::lattice.getStatus("current.focus.row", prefix = cprefix)
            column <- lattice:::lattice.getStatus("current.focus.column", prefix = cprefix)
		    switch(side,
              top = 
              {
			    if(row == 1 & column == 1)
                {				
			      #grid.text("  ", x = unit(-1, "npc"), y = unit(1.1, "npc"), just = "left") 
                  #grid.text("Legend", x = unit(-.1, "npc"), y = unit(1.1, "npc"), just = "right")
                  grid.text(c("recruitment","ssb"), x = c(.4, .6), hjust = c(1,0)) 
				}
				if (row == 1 & column > 1)
                  grid.text(c("", "Mesh\n","ICES areas\n")[column], 
				    x = unit(0.5, "npc"), y = unit(2, "npc"), hjust = 0.5, vjust = 0)
			  })
          })
		  
rec.p <-
levelplot(s.rec ~ x * y | factor(year), data = p.dat,
          colorkey = FALSE, as.table = TRUE, 
		  xlab = "", ylab = "", layout = c(3,2), strip = FALSE,
		  between = list(x = 0.1, y = 0.2),
		  par.settings = list(axis.line = list(col = "transparent")),
		  panel = panel.fun,
		  axis = function(side, ...) axis.fun(side, plot.label = c("a)", "Recruitment"), ...))

ssb.p <-
levelplot(s.ssb ~ x * y | factor(year), data = p.dat,
          colorkey = FALSE, as.table = TRUE, 
		  xlab = "", ylab = "", layout = c(3,2), strip = FALSE,
		  between = list(x = 0.1, y = 0.2),
		  par.settings = list(axis.line = list(col = "transparent")),
		  panel = panel.fun,
		  axis = function(side, ...) axis.fun(side, plot.label = c("b)", "SSB"), ...))
		  

cm.cols <- colorRampPalette(c("grey90","purple","cyan","yellow","red"))
cm.cols2 <- colorRampPalette(c("red", "pink","white"))
breaks <- c(seq(0, .4, length = 201), seq(.4, 1, length = 201)[-1]) 
cols <- c(cm.cols(200), cm.cols2(200))
		  
		  
		  
tiff(paste("figures/smooth.rec.ssb.tiff",sep=""), 8, 12, res = 500, units = "in")
  plot(legend.p, position = c(0.15, .57, 1, 1), more = TRUE)
  plot(rec.p, position = c(0.15, .37, 1, .8), more = TRUE)
  plot(ssb.p, position = c(0.15, 0, 1, .43), more = FALSE)	  
dev.off()

