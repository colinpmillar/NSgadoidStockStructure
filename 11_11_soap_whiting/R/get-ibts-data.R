###########################################################################
## Author: Colin Millar
## This code comes with no warranty or guarantee of any kind.
###########################################################################

setwd("~/work/11_11_soap_whiting")

########################
## fit function inputs: 
########################

## data preparation - a data frame directly from DATRAS with extra columns for observations (e.g. log recruitment)

ibts.fnames <- paste("data", c("NSibts.csv", "WCibts.csv"), sep="/")
ibts <- do.call(rbind, lapply(ibts.fnames, read.table, sep=",", header = TRUE, stringsAsFactors = FALSE))

# select only whiting
ibts <- subset(ibts, Species == "Merlangius merlangus" & Quarter == 1)

# get rectuits
ibts $ rec <- ibts $ Age_1

## get ssb
# fill a data frame with weight and maturity data
weight <- mat <- expand.grid(Survey = unique(ibts $ Survey), stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

weight[paste("Age_", 0:10, sep="")] <- rbind(c(0   ,0    ,0.179,0.239,0.285,0.310,0.348,NA  ,NA  ,NA  ,NA   ),  # NS Whi
					   					                       c(0   ,0    ,0.177,0.214,0.248,0.287,0.299,NA  ,NA  ,NA  ,NA   ))  # WC Whi

mat[paste("Age_", 0:10, sep="")] <- rbind(c(0  ,0   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ),  # NS Whi
										                      c(0  ,0   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ))  # WC Whi
										  
ibts $ ID <- 1:nrow(ibts)
weight <- merge(ibts[c("Survey","ID")], weight, all = TRUE)
weight <- weight[order(weight $ ID),] 
mat <- merge(ibts[c("Survey","ID")], mat, all = TRUE)
mat <- mat[order(mat $ ID),] 
					
cnames <- paste("Age_",0:10,sep="")
					
ibts $ ssb <- rowSums(weight[cnames] * mat[cnames] * ibts[cnames], na.rm = TRUE)
							
# my personal preference
names(ibts) <- tolower(names(ibts))

#ibts <- ibts[c("year", "survey", "shootlat", "shootlon", "rec", "ssb")]

save(ibts, file="run/whiting-ibts.RData")
