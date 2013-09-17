
################
## Utility function
################

require(INLA)

create_surface <-
function(limits, max.edge = 1)
{
  segm_bndry = inla.mesh.segment(loc = limits[[1]])
  segm_inners = lapply(limits[-1], function(x) inla.mesh.segment(loc = x[nrow(x):1,]))
  
  inla.mesh.create(boundary = c(list(segm_bndry), segm_inners), refine = list(max.edge = max.edge))
}

cm.cols <- colorRampPalette(c("grey90","purple","cyan","yellow","red"))
cm.cols2 <- 
function(i) 
{
  cols1 <- colorRampPalette(c("grey90","purple","cyan","yellow","red"))
  cols2 <- colorRampPalette(c("red", "pink","white"))
  c(cols1(floor(i/2)), cols2(ceiling(i/2)))
}

lonlat3D <-
function(xy)
{
  lon <- xy[,1]
  lat <- xy[,2]
  cbind(sin((lon/180)*pi)*cos((lat/180)*pi),
		sin((lat/180)*pi),
		cos((lon/180)*pi)*cos((lat/180)*pi))
}

lonlat2D <- # works okay for the quadrant of the globe Europe is in
function(xyz)
{
  x <- xyz[,1]
  y <- xyz[,2]
  z <- xyz[,3]
  cbind(asin(x / cos(asin(y))),
		asin(y)  ) * 180 / pi
}


################
## The fit function
################


fit_surface <-
function(Data, mesh, init.mode = NULL, verbose = FALSE)
{

  #Convenient definitions - number of vertices and number of observations

  nV = mesh $ n
  # Make the spde - simple non-intrinsic Matern
  spde = inla.spde.create(mesh, model = "matern")
  frm = y ~ f(field, model = spde) - 1

  # trim data
  Data_proj <- inla.mesh.project(mesh, Data $ loc)
  Data <- list(loc = Data $ loc[Data_proj $ ok,], y = Data $ y[Data_proj $ ok])

  # create a vector of data + NAs (one for each node) for prediction and its design matrix
  y <- c(Data $ y, rep(NA, nV)) 
  A = rBind(inla.mesh.project(mesh, Data $ loc) $ A, Matrix:::sparseMatrix(i = 1:nV, j = 1:nV, x = rep(1,nV)) )


 if (!is.null(init.mode)) {

  result = inla(frm, data = list(y = y, field = 1:nV), family = "normal", 
          control.predictor = list(A = A, compute = TRUE), num.threads = 4, 
          control.mode = list(theta = init.mode, restart=TRUE),
          verbose = verbose)

  } else {

  result = inla(frm, data = list(y = y, field = 1:nV), family = "normal",  
                control.predictor=list(A = A, compute=TRUE),
	    		num.threads = 1, verbose = verbose)
  }
  
  result $ summary.linear.predictor[length(Data $ y) + 1:nV,]
}

################
## plotting functions
################

plot.data <-
function(loc, mesh, limits, land, depths, grd, fname = NULL)
{
  tofile <- !is.null(fname)
  if (tofile) tiff(paste(fname, ".tiff", sep=""), 8, 8, res = 500, units = "in")

  plot.new()
  if (!missing(mesh)) 
    plot.window(xlim = range(mesh$loc[, 1]), ylim = range(mesh$loc[, 2]), "")
  else 
    plot.window(xlim = range(loc[, 1]), ylim = range(loc[, 2]), "")
  if (!missing(land)) tmp <- lapply(land, polygon, col = "darkseagreen")
  if (!missing(depths)) tmp <- lapply(depths, lines, col = grey(.7))
  if (!missing(mesh)) plot(mesh, col = grey(0.4), add = TRUE)
  if (!missing(limits)) tmp <- lapply(limits, function(xy) polygon(xy, border = grey(0.3), lwd = 2))
  if (!missing(grd)) lines(grd, col = grey(0.4))
  points(loc, pch = 16, col = "red", cex = 0.5)
  
  if (tofile) dev.off()
}



plot.fit <-
function(z, proj, colfun = cm.cols, ncols = 50, ncontours = 4, 
         limits, land, breaks, clevels, zmax = max(z, na.rm = TRUE), fname = NULL)
{
  if (missing(breaks))
  {
    breaks <- c( quantile(z, c(0:(ncols-1)/ncols), na.rm = TRUE), zmax)
	  breaks <- round(signif(breaks, 2), 2)
  }
  if (missing(clevels)) 
  {
    clevels <- quantile(z, c(1:ncontours/(ncontours+1)), na.rm = TRUE)
	  clevels <- round(signif(clevels, 2), 2)
  }
  cols <- colfun(length(breaks) - 1)
  
  tofile <- !is.null(fname)
  if (tofile) tiff(paste(fname, ".tiff", sep=""), 8, 8, res = 500, units = "in")

  plot(0, 0, type = "n", ylab = "", xlab = "", las = 1, xlim = range(proj $ x), ylim = range(proj $ y))
  if (!missing(land)) tmp <- lapply(land, polygon, col = "darkseagreen", border = grey(0.8))
  image(proj $ x, proj $ y, z, col = cols, breaks = breaks, add = TRUE)
  if (ncontours > 0) contour(proj $ x, proj $ y, z, add = TRUE, levels = clevels, labels = clevels, labcex = .4)
  if (!missing(limits)) tmp <- lapply(limits, function(xy) polygon(xy, border = grey(0.3)))
  
  if (tofile) dev.off()
}


#########
## extra functions
##########


inla.mesh.lattice.cm <-
function (x = seq(0, 1, length.out = 2), y = seq(0, 1, length.out = 2), z = NULL, 
          dims = (inla.ifelse(is.matrix(x), dim(x), c(length(x), length(y)))), units = NULL) 
{
    units = match.arg(units, c("default", "longlat", "longsinlat"))
    if (missing(x) && !missing(dims)) {
        x = seq(0, 1, length.out = dims[1])
    }
    if (missing(y) && !missing(dims)) {
        y = seq(0, 1, length.out = dims[2])
    }
    dims = as.integer(dims)
    if (is.matrix(x)) {
        if (!identical(dims, dim(x)) || !identical(dims, dim(y)) || 
            (is.matrix(z) && !identical(dims, dim(z)))) 
            stop("The size of matrices 'x', 'y', and 'z' must match 'dims'.")
        loc = cbind(as.vector(x), as.vector(y), as.vector(z))
        x = NULL
        y = NULL
    }
    else {
        if (!identical(dims[1], length(x)) || !identical(dims[2], 
            length(y))) 
            stop(paste("The lengths of vectors 'x' and 'y' (", 
                length(x), ",", length(y), ") must match 'dims' (", 
                dims[1], ",", dims[2], ").", sep = ""))
        loc = (cbind(rep(x, times = dims[2]), rep(y, each = dims[1])))
    }
    if (!is.double(loc)) 
        storage.mode(loc) = "double"
    if (identical(units, "longlat")) {
        loc = lonlat3D(loc)
    }
    else if (identical(units, "longsinlat")) {
        coslat = sapply(loc[, 2], function(x) sqrt(max(0, 1 - 
            x^2)))
        loc = (cbind(cos(loc[, 1] * pi/180) * coslat, sin(loc[, 
            1] * pi/180) * coslat, loc[, 2]))
    }
    segm.idx = (c(1:(dims[1] - 1), dims[1] * (1:(dims[2] - 1)), 
        dims[1] * dims[2] - (0:(dims[1] - 2)), dims[1] * ((dims[2] - 
            1):1) + 1))
    segm.grp = (c(rep(1L, dims[1] - 1), rep(2L, dims[2] - 1), 
        rep(3L, dims[1] - 1), rep(4L, dims[2] - 1)))
    segm = (inla.mesh.segment(loc = loc[segm.idx, , drop = FALSE], 
        grp = segm.grp, is.bnd = TRUE))
    lattice = list(dims = dims, x = x, y = y, loc = loc, segm = segm)
    class(lattice) = "inla.mesh.lattice"
    return(lattice)
}

inla.mesh.projector.cm <-
function (mesh, loc = NULL, lattice = NULL, 
          xlim = range(mesh$loc[, 1]), ylim = range(mesh$loc[, 2]), 
          dims = c(100, 100), projection = NULL) 
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")
    if (missing(loc) || is.null(loc)) {
        if (missing(lattice) || is.null(lattice)) {
            if (identical(mesh$manifold, "R2")) {
                units = "default"
                x = seq(xlim[1], xlim[2], length.out = dims[1])
                y = seq(ylim[1], ylim[2], length.out = dims[2])
            }
            else if (identical(mesh$manifold, "S2")) {
                projection = match.arg(projection, c("longlat", 
                  "longsinlat"))
                units = projection
                if (missing(xlim) || is.null(xlim)) {
                  xlim = c(-180, 180)
                }
                if (missing(ylim) || is.null(ylim)) {
                  ylim = c(-90, 90)
                }
                x = seq(xlim[1], xlim[2], length.out = dims[1])
                if (identical(projection, "longlat")) {
                  y = seq(ylim[1], ylim[2], length.out = dims[2])
                }
                else {
                  y = (seq(sin(ylim[1] * pi/180), sin(ylim[2] * 
                    pi/180), length.out = dims[2]))
                }
            }
            lattice = (inla.mesh.lattice.cm(x = x, y = y, units = units))
        }
        else {
            dims = lattice$dims
            x = lattice$x
            y = lattice$y
        }
        proj = inla.mesh.project(mesh, lattice$loc)
        projector = list(x = x, y = y, lattice = lattice, loc = NULL, 
            proj = proj)
        class(projector) = "inla.mesh.projector"
    }
    else {
        proj = inla.mesh.project(mesh, loc)
        projector = list(x = NULL, y = NULL, lattice = NULL, 
            loc = loc, proj = proj)
        class(projector) = "inla.mesh.projector"
    }
    return(projector)
}


#########
## plotting functions
##########


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
		  
panel.fun <- function(x, y, z, subscripts, years, ...)
		  {
			cprefix <- lattice:::lattice.getStatus("current.prefix")
            row <- lattice:::lattice.getStatus("current.focus.row", prefix = cprefix)
            column <- lattice:::lattice.getStatus("current.focus.column", prefix = cprefix)
			tmp <- lapply(land, panel.polygon, col = "darkseagreen", border = FALSE)
			panel.levelplot(x, y, z, at = breaks, col.regions = cols, subscripts = subscripts)
			tmp <- lapply(limits, function(xy) panel.polygon(xy, border = grey(0.3)))
			#llines(split.line $ x, split.line $ y, col = "blue")
			grid.text(years[column + (row-1)*3], x = unit(0.025, "npc"), y = unit(.95, "npc"), just = "left")
		  }


panel.legend <- function(x, y, z, subscripts, ...)
{
  cprefix <- lattice:::lattice.getStatus("current.prefix")
  column <- lattice:::lattice.getStatus("current.focus.column", prefix = cprefix)
	if (column == 2) {
		 tmp <- lapply(land, panel.polygon, col = "darkseagreen", border = FALSE)
		 t.sub = 1:nrow(mesh $ graph $ tv)
		 idx = cbind(mesh $ graph $ tv[t.sub, c(1:3, 1), drop = FALSE], NA)
     x = mesh$loc[t(idx), 1]
     y = mesh$loc[t(idx), 2]
     llines(x, y, col = grey(0.5), lwd = 1)
		 tmp <- lapply(limits, function(xy) panel.polygon(xy, border = grey(0.3)))
	} else if (column == 3) {
	  lpolygon(ICES.areas $ x, ICES.areas $ y)
	  tmp <- lapply(land, lpolygon, col = "darkseagreen", border = FALSE)
		grid.text("IV", x = 2.72, y = 56.62, default.units = "native", just = "centre")
	} else if (column == 1) {
	  # do scale...
		panel.levelplot(scale.dat $ x, scale.dat $ y, scale.dat $ z, 
			              at = breaks, col.regions = cols, subscripts = subscripts)
		at <- seq(0, 1, by = 0.2)
		rat <- pretty(ifrec(at))
		sat <- pretty(ifssb(at))
		
		# dont plot too hight to avoid overlap with legend text
		rat <- rat[frec(rat) < 0.95]
		sat <- sat[fssb(sat) < 0.95]
		
		grid.text(sprintf(rat, fmt = "%5.0f"), 
			        x = rep(min(scale.dat $ x) - .3, length(rat)), 
			        y = ifz(frec(rat)), just = "right", default.units = "native", gp = gpar(cex = 0.7))
	  grid.text(sprintf(sat, fmt = "%4.0f"), 
			        x = rep(max(scale.dat $ x) + .3, length(sat)), 
			        y = ifz(fssb(sat)), hjust = 0, default.units = "native", gp = gpar(cex = 0.7))
	}
}

axis.legend <- function(side, ...)
{
  cprefix <- lattice:::lattice.getStatus("current.prefix")
  row <- lattice:::lattice.getStatus("current.focus.row", prefix = cprefix)
  column <- lattice:::lattice.getStatus("current.focus.column", prefix = cprefix)
	switch(side,
    top = {
		  if(row == 1 & column == 1) {				
			  #grid.text("  ", x = unit(-1, "npc"), y = unit(1.1, "npc"), just = "left") 
        #grid.text("Legend", x = unit(-.1, "npc"), y = unit(1.1, "npc"), just = "right")
        #grid.text(c("Recruitment (No. h-1)","SSB (kg h-1)"), x = c(.4, .6), hjust = c(1,0))
        grid.text(expression(paste("Recruitment (no. ", h^-1, ")"), paste("SSB (kg ", h^-1, ")")), x = c(.4, .6), hjust = c(1,0))
			}
			if (row == 1 & column > 1) {
        grid.text(c("", "Mesh\n","ICES areas\n")[column], 
				          x = unit(0.5, "npc"), y = unit(2, "npc"), hjust = 0.5, vjust = 0)
			}
	  })
}

