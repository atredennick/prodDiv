# This is an "R" script documenting work associated with the following
# publication: 
# Fraser, L.H., Pither, J., Jentsch, A., Sternberg, M., Zobel, M., and the 
# HerbDivNet # global network.  2015. Worldwide evidence of a unimodal 
# relationship between productivity and plant species richness.  Science.

# P.I.: Lauchlan Fraser (Thompson Rivers University) (lfraser@tru.ca)

# Script author: Jason Pither (University of British Columbia, Okanagan campus)
# (jason.pither@ubc.ca)

# ---------------------------------------------------------
# THIS SCRIPT CONDUCTS READS IN AND FORMATS THE DATA
# ---------------------------------------------------------

# ---------------------------------------------------------
# REQUIRED FILES:
# fraser_plotdata.csv
# ---------------------------------------------------------

# ---------------------------------------------------------
# read in data, and do initial data cleaning
# ---------------------------------------------------------

library(rdryad)
alldata <- dryad_getfile("http://datadryad.org/bitstream/handle/10255/dryad.90398/fraser_plotdata.csv?sequence=1")

# first 10 columns are summary data
alldata.summary <- alldata[,c("grid","original.name","community.type",
                              "country","pi","plot","biomass","litter",
                              "tot.bio","sr")]  

# column headings are self-explanatory; but note that "pi" is what we equate with "site" throughout the paper

# exclude grids that had not litter measurements conducted in all 64 quadrats
complete.grid.numbers <- setdiff(unique(alldata.summary$grid),
                                 names(table(alldata.summary[is.na(alldata.summary$litter),"grid"])
                                       [table(alldata.summary[is.na(alldata.summary$litter),"grid"]) == 64]))

num.complete.grids <- length(complete.grid.numbers)

# now we have 151 grids with complete or mostly complete measurements
# (10 have between 1 and 4 quadrats missing litter)

# now subset data to only those grids
alldata.litter <- alldata.summary[!is.na(match(alldata.summary$grid,
                                               as.numeric(complete.grid.numbers))),]

# there are some quadrats with zero species, and two (7023, 7030)
# in grid 110 that have valid litter biomass (zero) but no live species;
# these remain as VALID data points

# put "NA" in all other quadrats
na.rownames <- row.names(alldata.litter[alldata.litter$sr  ==  0,][setdiff(row.names(alldata.litter[alldata.litter$sr  ==  0,]),c("7023","7030")),])

# remaining quadrats with neither biomass nor species can be "NA"
alldata.litter[na.rownames,"biomass"] <- NA
alldata.litter[na.rownames,"litter"] <- NA
alldata.litter[na.rownames,"tot.bio"] <- NA
alldata.litter[na.rownames,"sr"] <- NA

# create a log10 total biomass (+1) column
alldata.litter$log10.tot.bio <- log10(alldata.litter$tot.bio+1) 

# CHECK which if any quadrats are missing live biomass:
alldata.litter[is.na(alldata.litter$tot.bio),]

## These quadrats must be excluded from analyses.

good.data <- alldata.litter[!is.na(alldata.litter$tot.bio),];

# there are 33 quadrats excluded, so the null DF for the global models is
# 9664 (num rows in alldata.litter)  minus 33 = 9631 (num rows in good.data)

# treat "grid" as a factor for analyses
good.data$grid <- factor(good.data$grid)

# same with "pi", which we consider equivalent to "site"
good.data$pi <- as.character(good.data$pi)
good.data$pi <- as.factor(good.data$pi)