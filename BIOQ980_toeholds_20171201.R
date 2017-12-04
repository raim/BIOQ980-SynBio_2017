library(platexpress)


## READING IN DATA: this is often the most time-consuming
## step, since softwares and users always change their export
## formats and no standard format for platereader data exists
## HERE: oocalc > export to csv; unticking "save data as shown" (usually
## recommmended) messes up times and introduces jumps! 

raw <- readPlateData(files="BIOQ980_toeholds_20171201_synergy.csv",
                     type="Synergy", skip=56, time.conversion=1/60,
                     time.format="%H:%M:%S",sep=",")

viewPlate(raw)

plate <- readPlateMap("2017_12_01_toeholds_plate_map.csv",sep=",",fsep=";",
                      fields=c("group","IPTG"), blank.id="LB")

g1 <- getGroups(plate,by="group")
g2 <- getGroups(plate,by=c("group","IPTG"))

viewGroups(raw,g1,g2)

raw <- prettyData(raw, dids=c(OD="600",GFP="GFP_50:480,520"))
data <- cutData(raw, rng=c(0,700))

viewGroups(data, g1,g2)


fl <- getData(data, "GFP")
od <- getData(data, "OD")
data <- addData(data, ID="FL/OD", dat=fl/od)

viewGroups(data, g1, g2, dids=c("OD","FL/OD"))

flod <- getData(data, "FL/OD")

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

gr1.mn <- apply(flod[,gr1.uni], 1, mean, na.rm=TRUE)
gr4.mn <- apply(flod[,gr4.uni], 1, mean, na.rm=TRUE)

uni <- flod
uni[,gr1.all] <- uni[,gr1.all]/gr1.mn
uni[,gr4.all] <- uni[,gr4.all]/gr4.mn

data <- addData(data, ID="FL/OD/unind.", dat=uni, col=wavelength2RGB(450),
                replace=TRUE)

viewGroups(data, g1,g2, dids=c("OD","FL/OD/unind."))

### OR: simple version by Louis

## interpolate to OD

od.data <- interpolatePlateData(data, "OD")

results <- boxData(od.data, did="FL/OD/unind.", rng=.6, groups=g2, type="bar")
