RD = R
DD = data
BD = run
FD = fig

.PHONY = clean data run plot show

clean:
	rm -f $(BD)/*
	rm -f $(FD)/*

data: $(BD)/whiting-ibts.RData $(BD)/boundary-data.RData
$(BD)/whiting-ibts.RData: $(DD)/NSibts.csv $(DD)/WCibts.csv $(RD)/get-ibts-data.R
	echo 'source("R/get-ibts-data.R")' | R --vanilla --slave
$(BD)/boundary-data.RData: $(DD)/contours.rdata $(DD)/land.rdata $(DD)/limits.rdata $(RD)/get-boundary-data.R
	echo 'source("R/get-boundary-data.R")' | R --vanilla --slave

run: $(BD)/output.RData
$(BD)/output.RData: $(BD)/whiting-ibts.RData $(BD)/boundary-data.RData $(RD)/run-model.R
	echo 'source("R/run-model.R")' | R --vanilla --slave

plot: $(BD)/output.RData $(RD)/paper-plot.R
	echo 'source("R/paper-plot.R")' | R --vanilla --slave

show: $(FD)/Rplot001.png
	eog fig/*.png &
  
