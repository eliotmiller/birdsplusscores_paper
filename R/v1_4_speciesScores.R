library(clootl)
library(phytools)
library(dplyr)
library(viridis)
library(Rphylopars)
library(ebirdst)
library(picante)

#set the working directory. this should be the directory that contains the data
#and R folders. if you got here by double-clicking the R script, the following
#will probably work
setwd("../")

##### there are toggles throughout script to go to an unlogged form
# of various variables if you want to experiment with that. set it back to
# the logged form when you are done

############################################################################
################### Assemble the vulnerability factors #####################
############################################################################

#try to obtain human tolerance data first. following the Marjakangas et al.
#2024 paper (https://onlinelibrary.wiley.com/doi/10.1111/geb.13816), I downloaded the dataset
#from here https://datadryad.org/stash/dataset/doi:10.5061/dryad.83bk3jb08 on 18 Mar 2024
#These authors seem to prefer conservative HTI, so use that column. the file here is simply
#a CSV of the data sheet from the downloaded XLSX file
htiDat <- read.csv("data/bird_species_tolerances_dataset_240228.csv")

#pull a tree. note that this is 2023 taxonomy, and the HTI file is 2022 taxonomy.
tree <- extractTree("all_species", label_type="scientific", taxonomy_year=2023)

#make the tree ultrametric. this is code I've carried along through various 
#scripts from Jonathan Chang.
tree <- reorder(tree, "postorder")
e1 <- tree$edge[, 1]
e2 <- tree$edge[, 2]
EL <- tree$edge.length
N <- Ntip(tree)
ages <- numeric(N + tree$Nnode)
for (ii in seq_along(EL)) {
  if (ages[e1[ii]] == 0) {
    ages[e1[ii]] <- ages[e2[ii]] + EL[ii]
  }
  else {
    recorded_age <- ages[e1[ii]]
    new_age <- ages[e2[ii]] + EL[ii]
    if (recorded_age != new_age) {
      EL[ii] <- recorded_age - ages[e2[ii]]
    }
  }
}

tree$edge.length <- EL
tree <- ladderize(tree)

#add underscores and extract any species that match 2023 taxonomy. later, bring in Avibase
#taxon codes and properly update this file to 2023 taxonomy, then do all analyses in that
#taxonomy version
htiDat$underscores <- sub(" ", "_", htiDat$Scientific_name_eBird)

#there are continent-specific HTI scores for species on multiple continents. take the
#median of these for now. realistically might want to exclude continents with low
#abundance in the future
forSummary <- htiDat[htiDat$underscores %in% tree$tip.label, c("underscores","Conservative_HTI")]
forSummary <- group_by(forSummary, underscores)
htiFinal <- as.data.frame(summarize(forSummary, hti=median(Conservative_HTI)))

#interestingly, did not lose any species there. looks like they already back-calculated it to 2021
#plop this into a contMap and see how it looks. for now just plot on a pruned tree. impute
#missing values later.
x <- htiFinal$hti
names(x) <- htiFinal$underscores

pruned <- drop.tip(tree, setdiff(tree$tip.label, names(x)))

# obj <- contMap(pruned, x, fsize=0.05, outline=FALSE, lwd=0.5, res=200, plot=FALSE)
# n <- length(obj$cols)
# obj$cols[1:n] <- plasma(n)
# pdf(file="outputs/hti_30May2024.pdf", height=80, width=40)
# plot(obj, fsize=0.05, lwd=0.5, outline=FALSE)
# dev.off()

#now try bringing in the range size data. this is an older file from Uri Roll
#that was used in some publications by those authors
rsDat <- read.csv("data/birds_gbd_bioclim.csv")

#drop to just the two columns of interest
rsDat <- rsDat[,c("ScinameGBD","Bio1_count")]

#and bring in some taxonomic translation tables from AVONET to help convert.
#the column Species3 here is the Jetz tree, Species1 is BirdLife, Species2
#is eBird
BLtoBT <- read.csv("data/BirdLife-BirdTree-crosswalk.csv")
BLtoeBird <- read.csv("data/BirdLife-eBird-crosswalk.csv")

#rsDat isn't clearly in either taxonomy, but it matches better to BT (bird tree), so use that to
#convert to BL, and then to eBird
merge1 <- merge(rsDat, BLtoBT[,c("Species1","Species3")], by.x="ScinameGBD", by.y="Species3")
merge1 <- merge1[!duplicated(merge1$Species1),]
merge2 <- merge(merge1, BLtoeBird[,c("Species1","Species2")], by.x="Species1", by.y="Species1")
merge2 <- merge2[!duplicated(merge2$Species2),]

#create an underscore column
merge2$underscores <- sub(" ", "_", merge2$Species2)

#there are a few non-species-level taxa here. drop these
merge2 <- merge2[merge2$underscores %in% tree$tip.label,]

#set this result aside
rangeSize <- data.frame(underscores=merge2$underscores, bl.range.size=merge2$Bio1_count)

#also bring in migration distance. comes from DOI: 10.1111/jbi.13700
migDist <- read.csv("data/jbi13700-sup-0002-datas1.csv")

#there is a weird extra row here that screws up merges later
migDist <- migDist[migDist$Sp.Scien.jetz != "",]

#try to get taxonomy matched
migDist$no.underscores <- gsub("_", " ", migDist$Sp.Scien.jetz)

#match jetz to birdlife, then birdlife to ebird
merge3 <- merge(migDist, BLtoBT[,c("Species1","Species3")], by.x="no.underscores", by.y="Species3")
merge3 <- merge3[!duplicated(merge3$Species1),]
merge4 <- merge(merge3, BLtoeBird[,c("Species1","Species2")], by.x="Species1", by.y="Species1")
merge4 <- merge4[!duplicated(merge4$Species2),]

#create an underscore column
merge4$underscores <- sub(" ", "_", merge4$Species2)

#there are a few non-species-level taxa here. drop these
merge4 <- merge4[merge4$underscores %in% tree$tip.label,]

#there are also some NAs. drop those too. we are using the "ALL" migration distance
#column here
merge4 <- merge4[!is.na(merge4$distance_quanti_ALL),]

#set these results aside
migResults <- data.frame(underscores=sub(" ", "_", merge4$Species2),
                         mig.distance=merge4$distance_quanti_ALL)

#load in the eBird migration distance. you calculate this yourself using the 
#script ebirdMigrationDistance.R. this is 2023 taxonomy
eMig <- read.csv("data/eBirdMigDist_1Dec2024.csv")

# load the taxonomy
tax <- read.csv("data/ebird_taxonomy_v2023.csv")

# merge and add an underscores column
eMig <- merge(eMig, tax[,c("SPECIES_CODE","SCI_NAME")], by.x="species.code", by.y="SPECIES_CODE")
eMig$underscores <- sub(" ", "_", eMig$SCI_NAME)

# make a nicer name for the data column
names(eMig)[2] <- "ebird.mig.dist"

# last bring in eBird S&T. you pre-processed these products in a separate script
ebird <- read.csv("data/processedST.csv")

# merge with taxonomy and cut to species
ebird <- merge(ebird, tax[,c("SPECIES_CODE","SCI_NAME", "CATEGORY")], by.x="species",
               by.y="SPECIES_CODE")
ebird <- ebird[ebird$CATEGORY=="species",]

# give underscores
ebird$underscores <- sub(" ", "_", ebird$SCI_NAME)

############################################################################
################### Phylogenetic and functional uniqueness #################
############################################################################

#bring in avonet
avo <- read.csv("data/AVONET-ebirdTax.csv")

#set the variables aside that you want, then transform as necessary
avo <- avo[,c("Species2","Beak.Length_Culmen","Beak.Length_Nares","Beak.Width","Beak.Depth",
              "Tarsus.Length","Wing.Length","Secondary1","Tail.Length","Mass")]

# log these
temp <- apply(avo[,2:dim(avo)[2]], 2, log)
temp2 <- avo[,2:dim(avo)[2]]
avo <- data.frame(ebird.old=avo$Species2, temp)

#run the PCA
pca <- prcomp(avo[,2:dim(avo)[2]], center = TRUE, scale. = TRUE)

#find the centroid, then calculate the distance of each species from that centroid
centroid <- colMeans(pca$x)
dists <- mahalanobis(x=pca$x, center=centroid, cov=cov(pca$x))

#add names and check the results out. results look good. now plot results
names(dists) <- avo$ebird.old

distResults <- data.frame(underscores=sub(" ", "_", avo$ebird.old), avo, sp.dist=dists)
distResults <- distResults[distResults$underscores %in% tree$tip.label,]
distResults$ebird.old <- NULL

# calculate these values. takes < 30s
system.time(edgeScores <- evol.distinct(tree, "equal.splits"))

# get names right for later
names(edgeScores) <- c("underscores", "edge.score")

# load the processed geographic information as well. the script & data to create this
# file is included as well
ebirdPhyloMorpho <- read.csv("data/ebirdPhyloMorpho_3Mar2025.csv")
names(ebirdPhyloMorpho) <- c("name","gen.field","morph.field")
tax$underscores <- sub(" ", "_", tax$SCI_NAME)
ebirdPhyloMorpho <- merge(ebirdPhyloMorpho,
                          tax[,c("SPECIES_CODE","underscores")],
                          by.x="name", by.y="SPECIES_CODE")

############################################################################
######################### Conservation status ##############################
############################################################################

# now load in the IUCN threat data from J Gerbracht
iucn <- read.csv("data/Clements2023-IUCN.csv")

# convert these threat categories to ordinal scale. your current conversion table
# is here
ordConv <- read.csv("data/threats_to_ordinal_v1-3.csv")

# change a name for ease later
names(ordConv)[2] <- "iucn.ordinal"

iucn <- merge(iucn, ordConv, by.x="value", by.y="category")

# add an underscores column
iucn$underscores <- sub(" ", "_", iucn$sci_name)

# load the recovered IUCN codes too. this was some manual taxonomy matching i did
newIUCN <- read.csv("data/recoveredIUCN.csv")
newConv <- read.csv("data/qualitativeIUCN_v1-3.csv")
newIUCN$underscores <- sub(" ", "_", newIUCN$scientificName)
newIUCN <- newIUCN[,c("underscores","redlistCategory")]

# merge in the new ordinal scores
newIUCN <- merge(newIUCN, newConv, by.x="redlistCategory", by.y="qualitative")
names(newIUCN)[3] <- "iucn.ordinal"

# bind these up
iucn <- iucn[,c("underscores","iucn.ordinal")]
iucn <- rbind(iucn, newIUCN[,c("underscores","iucn.ordinal")])

# load the ACAD scores
acadRaw <- read.csv("data/ACAD Global 2024.05.23.csv")

# lots of extra rows, delete
acadRaw <- acadRaw[!is.na(acadRaw$CCS.max),]

# grab what you need
acad <- data.frame(code2024=acadRaw$X2024.Species.Code,
                   sci.name=acadRaw$Scientific.Name, max.acad=acadRaw$CCS.max)

# look at what doesn't match (this is 2024 taxonomy)
issues <- setdiff(acad$code2024, tax$SPECIES_CODE)
acad[acad$code2024 %in% issues,]

# see if you can rescue these. we get the vast majority, just issues w/ western flycatcher
acad$code2023 <- acad$code2024

for(i in 1:dim(acad)[1])
{
  if(acad$code2024[i]=="#N/A")
  {
    # look this up and store match if there is one
    store <- tax$SPECIES_CODE[tax$SCI_NAME==acad$sci.name[i]]
    if(length(store) == 1)
    {
      acad$code2023[i] <- store
    }
    else
    {
      next()
    }
  }
}

# merge in underscores
acad <- merge(acad, tax[,c("SPECIES_CODE", "underscores")],
              by.x="code2023", by.y="SPECIES_CODE")

############################################################################
######################### Stack all the data into ##########################
########################### a single aligned DF ############################
############################################################################

unified <- merge(htiFinal, rangeSize, all=TRUE)
unified <- merge(unified, migResults, all=TRUE)
unified <- merge(unified, eMig[,c("underscores","ebird.mig.dist")], by="underscores", all=TRUE)
unified <- merge(unified, distResults, all=TRUE)
unified <- merge(unified, iucn[,c("underscores","iucn.ordinal")], all=TRUE)
unified <- merge(unified, ebird[,c("underscores","median.abund","range.size","trend")], all=TRUE)
unified <- merge(unified, edgeScores, all=TRUE)
unified <- merge(unified, ebirdPhyloMorpho[,c("gen.field","morph.field","underscores")], all=TRUE)
unified <- merge(unified, acad[,c("underscores","max.acad")], all=TRUE)

# merge in taxonomy and cut to just species
unified <- merge(unified,
                 tax[,c("underscores","PRIMARY_COM_NAME","SCI_NAME","CATEGORY","SPECIES_CODE")],
                 by="underscores")
unified <- unified[unified$CATEGORY=="species",]

# we aren't missing any species at this point, so no need to add rows for them.
# impute everything that's still missing, but transform as necessary first.
toImpute <- unified
toImpute$bl.range.size <- log(toImpute$bl.range.size + 1)
toImpute$range.size <- log(toImpute$range.size)
toImpute$mig.distance <- log(toImpute$mig.distance + 1)
toImpute$ebird.mig.dist <- log(toImpute$ebird.mig.dist + 1)
toImpute$sp.dist <- log(toImpute$sp.dist)
toImpute$median.abund <- log(toImpute$median.abund)
toImpute$edge.score <- log(toImpute$edge.score)
toImpute$edge.score <- log(toImpute$morph.field)

# cut to just names you want to pass in
dropNames <- c("underscores","PRIMARY_COM_NAME","SCI_NAME","CATEGORY","SPECIES_CODE")
tempNames <- toImpute$underscores
toImpute <- toImpute[,!(names(toImpute) %in% dropNames)]
toImpute <- data.frame(species=tempNames, toImpute)

# if you want to figure out what quantile a given value is in, you can do
# something like this
migECDF <- ecdf(toImpute$ebird.mig.dist)
migECDF(toImpute[toImpute$species=="Setophaga_cerulea","ebird.mig.dist"])

# save out what was imputed and what wasn't. use this to generate a "confidence"
# matrix, where species * trait information that was impute gets 10% of the 
# its normal weight in the averaging at the end
confDF <- toImpute
confDF[!is.na(confDF)] <- 1
confDF[is.na(confDF)] <- 0.05 #in version 1.1 you used 0.1 here
confDF$species <- toImpute$species

# run the imputation. OU gives best results. takes 330s
system.time(imputed <- phylopars(trait_data=toImpute, tree=tree, model="OU"))

#pull these out and further align data for BirdsPlus Index valuations
tempDF <- data.frame(imputed$anc_recon[1:length(tree$tip.label),])

min_max_scale <- function(x)
{
  (x - min(x)) / (max(x) - min(x))
}

#set up a finalDF
finalDF <- data.frame(species=row.names(tempDF))

# scale the variables from 0 to 1. flip range size so bigger values equal
# larger species scores. the same is also true of hti
finalDF$iucn.ordinal <- min_max_scale(tempDF$iucn.ordinal)
finalDF$max.acad <- min_max_scale(tempDF$max.acad)
finalDF$trend <- 1-min_max_scale(tempDF$trend)

finalDF$bl.range.size <- 1-min_max_scale(tempDF$bl.range.size)
finalDF$range.size <- 1-min_max_scale(tempDF$range.size)
finalDF$median.abund <- 1-min_max_scale(tempDF$median.abund)
finalDF$hti <- 1-min_max_scale(tempDF$hti)
finalDF$mig.distance <- min_max_scale(tempDF$mig.distance)
finalDF$ebird.mig.dist <- min_max_scale(tempDF$ebird.mig.dist)

finalDF$edge.score <- min_max_scale(tempDF$edge.score)
finalDF$sp.dist <- min_max_scale(tempDF$sp.dist)
finalDF$gen.field <- min_max_scale(tempDF$gen.field)
finalDF$morph.field <- min_max_scale(tempDF$morph.field)

# come up with a weights matrix for the averages by category
weightsDF <- finalDF

weightsDF$iucn.ordinal <- 1
weightsDF$max.acad <- 2
weightsDF$trend <- 3

weightsDF$bl.range.size <- 1
weightsDF$range.size <- 3
weightsDF$median.abund <- 4
weightsDF$hti <- 4
weightsDF$mig.distance <- 1
weightsDF$ebird.mig.dist <- 3

weightsDF$edge.score <- 1
weightsDF$sp.dist <- 1
weightsDF$gen.field <- 1
weightsDF$morph.field <- 1

# subset confidence DF to same columns and order as weights DF
confDF <- confDF[,names(confDF) %in% names(weightsDF)]
confDF <- confDF[,names(weightsDF)]

# get the rows in the same order as well!
row.names(confDF) <- confDF$species
row.names(weightsDF) <- weightsDF$species
confDF <- confDF[row.names(weightsDF),]

# multiply these matrices
finalWeights <- data.frame(species=confDF$species,
                           confDF[,2:dim(confDF)[2]] * weightsDF[,2:dim(weightsDF)[2]])
row.names(finalWeights) <- NULL

# now actually take the weighted averages and rescale again. also bring along the
# key variables so you can look at those in the future (but they don't contribute
# directly to the scores)
scoresDF <- finalDF[,names(finalWeights)]

conStatus <- c() 

for(i in 1:dim(scoresDF)[1])
{
  conStatus[i] <- weighted.mean(x=scoresDF[i,c("iucn.ordinal","max.acad","trend")],
                                w=finalWeights[i,c("iucn.ordinal","max.acad","trend")])
}

scoresDF$conservation.status <- min_max_scale(conStatus)

vulFactors <- c()

for(i in 1:dim(scoresDF)[1])
{
  vulFactors[i] <- weighted.mean(x=scoresDF[i,c("bl.range.size","range.size","median.abund","hti","mig.distance","ebird.mig.dist")],
                                w=finalWeights[i,c("bl.range.size","range.size","median.abund","hti","mig.distance","ebird.mig.dist")])
}

scoresDF$vulnerability.factors <- min_max_scale(vulFactors)

uni <- c()

for(i in 1:dim(scoresDF)[1])
{
  uni[i] <- weighted.mean(x=scoresDF[i,c("edge.score","sp.dist","gen.field","morph.field")],
                                 w=finalWeights[i,c("edge.score","sp.dist","gen.field","morph.field")])
}

scoresDF$uniqueness <- min_max_scale(uni)

scoresDF$bpi.score <- scoresDF$conservation.status +
  scoresDF$vulnerability.factors + scoresDF$uniqueness

#save this out
scoresDF <- merge(tax[,c("underscores","PRIMARY_COM_NAME","SPECIES_CODE","TAXON_ORDER")], scoresDF,
                  by.x="underscores", by.y="species")
scoresDF <- scoresDF[order(scoresDF$TAXON_ORDER),]
write.csv(scoresDF, "data/indexScores_v1-4.csv", row.names=FALSE)
write.csv(toImpute, "data/rawValues_v1-4.csv", row.names=FALSE)
write.csv(finalWeights, "data/finalWeights_v1-4.csv", row.names=FALSE)

