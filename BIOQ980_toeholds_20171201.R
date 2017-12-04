library(platexpress)


## READING IN DATA: this is often the most time-consuming
## step, since softwares and users always change their export
## formats and no standard format for platereader data exists
## HERE: oocalc > export to csv; unticking "save data as shown" (usually
## recommmended) messes up times and introduces jumps! 

## NOTE: that this function interpolates all data to the same
## time-points!
raw <- readPlateData(files="BIOQ980_toeholds_20171201_synergy.csv",
                     type="Synergy", skip=56, time.conversion=1/60,
                     time.format="%H:%M:%S",sep=",")

## inspect raw data!
viewPlate(raw)


## read plate layout
plate <- readPlateMap("2017_12_01_toeholds_plate_map.csv",sep=",",fsep=";",
                      fields=c("group","IPTG"), blank.id="LB")

## get groupings
g1 <- getGroups(plate,by="group")
g2 <- getGroups(plate,by=c("group","IPTG"))

## inspect grouped data
viewGroups(raw,g1,g2)

## cut to data range of interest
raw <- cutData(raw, rng=c(0,650))

## filter data, assign better names and colors
raw <- prettyData(raw, dids=c(OD="600",GFP="GFP_50:480,520"),
                  col=c("#000000",wavelength2RGB(520)))

## correct blank values (raw LB medium)
data <- correctBlanks(data=raw, plate=plate)

## inspect grouped data
viewGroups(data, g1,g2)


## Fluorescence per OD ~ proteins/cell
fl <- getData(data, "GFP")
od <- getData(data, "OD")
data <- addData(data, ID="FL/OD", dat=fl/od)

## DISCUSS: why does FL/OD (proteins per cell)
## vary over the growth phases!?
viewGroups(data, g1, g2, dids=c("OD","FL/OD"))

## FOLD-CHANGE induced over uninduced
flod <- getData(data, "FL/OD")

## get uninduced for G1 and G4 separately
## 1)  the hard way
gr1 <- plate[,"group"]=="G1"
gr4 <- plate[,"group"]=="G4"
gr1[is.na(gr1)] <- FALSE
gr4[is.na(gr4)] <- FALSE

unind <- plate[,"IPTG"]=="00"
unind[is.na(unind)] <- FALSE

gr1.all <- as.character(plate[gr1 , 1])
gr1.uni <- as.character(plate[gr1 & unind , 1])

gr4.all <- as.character(plate[gr4 , 1])
gr4.uni <- as.character(plate[gr4 & unind , 1])

## calculate mean for each time-point
gr1.mn <- apply(flod[,gr1.uni], 1, mean)
gr4.mn <- apply(flod[,gr4.uni], 1, mean)

## normalize separately for each group
uni <- flod
uni[,gr1.all] <- uni[,gr1.all]/gr1.mn
uni[,gr4.all] <- uni[,gr4.all]/gr4.mn

data <- addData(data, ID="FL/OD/unind.", dat=uni, col=wavelength2RGB(450),
                replace=TRUE)

## 2) the simple way:
## uninduced for each group are present already in grouping!
gr1.all <- g1[["G1"]] # NOTE the double [[]] for lists!
gr1.uni <- g2[["G1_00"]] 
gr4.all <- g1[["G4"]] 
gr4.uni <- g2[["G4_00"]] 
## normalize separately for each group
uni <- flod
uni[,gr1.all] <- uni[,gr1.all]/gr1.mn
uni[,gr4.all] <- uni[,gr4.all]/gr4.mn

data <- addData(data, ID="FL/OD/unind.", dat=uni, col=wavelength2RGB(450),
                replace=TRUE)


viewGroups(data, g1,g2, dids=c("OD","FL/OD/unind."))

## GET FINAL VALUE: fold-change at certain OD
## interpolate to OD - cut data to tighter range to avoid dual FL/OD values

od.data <- interpolatePlateData(cutData(data, rng=c(0,500)), "OD")
viewGroups(od.data, g1,g2, dids=c("FL/OD/unind."))

results <- boxData(od.data, did="FL/OD/unind.", rng=.8, groups=g2, type="bar")

## TODO: compare G1 vs. G4; 

## TODO: instead of at given OD, get the maximal FL/OD/unind. values,
## and record at which time/OD this is reached

## TODO: fit growth rates for each experiment/well
## by grofit and plot growth rate and maximal capacity vs. IPTG induction

## cut out first growth phase
gdat  <- cutData(data,rng=c(100,250))  
viewGroups(gdat, g1, g2, dids=c("OD","FL/OD"))


grodat <- data2grofit(gdat, did="OD", plate=plate, wells=unlist(g1))


