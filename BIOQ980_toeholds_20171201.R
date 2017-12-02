library(platexpress)


## READING IN DATA: this is often the most time-consuming
## step, since softwares and users always change their export
## formats and no standard format for platereader data exists
## HERE: oocalc > export to csv; unticking "save data as shown" (usually
## recommmended) messes up times and introduces jumps! 

raw <- readPlateData(files="BIOQ980_toeholds_20171201.csv",type="Synergy", skip=56, time.conversion=1/60, time.format="%H:%M:%S",sep=",")

viewPlate(raw)
