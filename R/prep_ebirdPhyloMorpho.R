library(picante)
library(clootl)
library(data.table)
library(foreach)
library(doParallel)
library(dplyr)
library(purrr)
library(Rphylopars)
library(phytools)
library(viridis)

#set the working directory. this should be the directory that contains the data
#and R folders. if you got here by double-clicking the R script, the following
#will probably work
setwd("../")

# now also calculate the phylogenetic uniqueness of every species in the global sense
tree <- extractTree("all_species", label_type="code", taxonomy_year=2023)

# make the tree ultrametric
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

# check how much memory it would require to make a giant phylo and morpho distance
# matrix, and whether you could parallelize across those
system.time(genDists <- cophenetic(tree))

# set diag to NA
diag(genDists) <- NA

# bring in avonet. this is in 2022 taxonomy
avo <- read.csv("data/AVONET-ebirdTax.csv")

# convert to 2023 taxonomy
tax2022 <- read.csv("data/2022.txt", sep="\t")
tax2022$avibase.id <- paste("AVIBASE", tax2022$bio_concept_code, sep="-")
avo <- merge(avo, tax2022[,c("species_code","avibase.id")],
             by.x="Avibase.ID2", "avibase.id", all.x=TRUE)
tax2023 <- read.csv("data/2023.txt", sep="\t")
tax2023$avibase.id <- paste("AVIBASE", tax2023$bio_concept_code, sep="-")
toMerge <- data.frame(new.code=tax2023$species_code,
                      avibase.id=tax2023$avibase.id, category=tax2023$category)
# there is some funky business here
toMerge <- toMerge[!duplicated(toMerge$avibase.id),]
toMerge <- toMerge[toMerge$category=="species",]

# continue with the merges
avo <- merge(avo, toMerge, by.x="Avibase.ID2", by.y="avibase.id", all.x=TRUE)

# drop information for anything that isn't a species now
avo <- avo[!is.na(avo$new.code),]

# set the variables aside that you want, then transform as necessary
newCode <- avo$new.code
avo <- avo[,c("Species2","Beak.Length_Culmen","Beak.Length_Nares","Beak.Width","Beak.Depth",
              "Tarsus.Length","Wing.Length","Secondary1","Tail.Length","Mass")]

# transform
temp <- apply(avo[,2:dim(avo)[2]], 2, log)
avo <- data.frame(code=newCode, temp)

# impute missing values
toImpute <- data.frame(species=avo$code, avo[,2:dim(avo)[2]])
missing <- setdiff(tree$tip.label, toImpute$code)
toBind <- data.frame(missing, matrix(nrow=length(missing), ncol=dim(toImpute)[2]-1))
names(toBind) <- names(toImpute)
toImpute <- rbind(toImpute, toBind)

# OU took 67s
system.time(imputed <- phylopars(trait_data=toImpute, tree=tree, model="OU"))

# pull the imputed data out
avoFinal <- data.frame(species=row.names(imputed$anc_recon)[1:length(tree$tip.label)],
                       imputed$anc_recon[1:length(tree$tip.label),])

# run the PCA. technically don't need PCA, but legacy code and it's fast
pca <- prcomp(avoFinal[,2:dim(avoFinal)[2]], center = TRUE, scale. = TRUE)

# calculate the distances between all species
toDist <- pca$x
system.time(morphDists <- as.matrix(dist(toDist, diag=TRUE, upper=TRUE)))

# set diag to NA
diag(morphDists) <- NA

# RAM and compute issues very low so far, let's try to parallelize over these.
# read in all of your eBird parsed data. this is in 2023 taxonomy
# identify all the files/species you will loop over
allFiles <- list.files("data/queried&processed/")
spp <- unlist(lapply(strsplit(allFiles, "_"), "[", 1))

#append the correct file prefix for loading quickly
allFiles <- paste("data/queried&processed/", 
                  allFiles, sep="")

# that worked. now try a parallel for loop
registerDoParallel(6)

# took 25s to read in
system.time(results <- foreach(theFile = 1:length(allFiles)) %dopar%
              {
                #read the file
                tempFile <- fread(allFiles[theFile])
                
                #if there are now valid rows, set to NA
                if(dim(tempFile)[1]==0)
                {
                  tempResult <- NA
                }
                
                else
                {
                  # else set it in place
                  tempResult <- as.data.frame(tempFile)
                }
                
                #return the result
                tempResult
              }
)

# now go through each species, and concatenate indices, which are index-year-month
# combinations, into just index-month. doing this for migratory reasons, previously
# was only retaining index to maximize data. not sure what impact will be of
# splitting data more finely now.

# now shifting back to drop both the month and year to speed up analysis
system.time(concat <- foreach(i = 1:length(results)) %dopar%
{
  if(is.null(dim(results[[i]])))
  {
    toStore <- NA
  }
  else
  {
    # come up with a new index
    temp <- strsplit(results[[i]]$index, "-")
    first <- unlist(lapply(temp, "[", 1))
    #second <- unlist(lapply(temp, "[", 3))
    #results[[i]]$new.index <- paste(first, second, sep="-")
    results[[i]]$new.index <- first
    
    # group and summarize (summing abundances)
    grouped <- group_by(results[[i]], new.index)
    toStore <- as.data.frame(summarize(grouped, n=sum(n)))
  }
  toStore
})

# remove results
rm(results)

# give names to concat
names(concat) <- spp

# drop NAs
concat <- concat[!is.na(concat)]

# combine the data frames with a new column for the species name. this comes from ChatGPT
combinedDF <- imap_dfr(concat, ~ mutate(.x, name = .y))

# identify all unique spp still remaining
spp <- unique(combinedDF$name)

# run this in parallel
registerDoParallel(8)

# loop over each species. identify all the unique indices it occurs in
system.time(spResults <- foreach(i = 1:length(spp)) %dopar%
{
  write.csv(data.frame(iteration=i), "data/progress.txt", row.names=FALSE)
  
  # identify all the indices it occurs in (don't think the unique is needed here
  # but just in case)
  indices <- unique(combinedDF$new.index[combinedDF$name==spp[i]])
  
  # set up index-level results vectors
  indexGenResult <- c()
  indexMorphResult <- c()
  spWeight <- c()
  
  # figure out which species are in the index and their abundances 
  # [could add some kind of statement
  # here to exclude indices with low sampling effort if you want]
  for(j in 1:length(indices))
  {
    otherSpp <- combinedDF[combinedDF$new.index==indices[j],]
    
    # how about if less than three species were recorded in that index you exclude
    if(dim(otherSpp)[1] < 3)
    {
      indexGenResult[j] <- NA
      indexMorphResult[j] <- NA
      spWeight[j] <- NA
    }
    else
    {
      # find the weighted distances between spp[i] and all other species. the
      # abundance of spp[i] doesn't matter here. do this in a loop to ensure
      # you retain the correct abundances of each otherSpp
      theWeight <- c()
      theGenDist <- c()
      theMorphDist <- c()
      
      for(k in 1:dim(otherSpp)[1])
      {
        # grab the distances and the weight
        theGenDist[k] <- genDists[spp[i], otherSpp$name[k]]
        theMorphDist[k] <- morphDists[spp[i], otherSpp$name[k]]
        theWeight[k] <- otherSpp$n[k]
      }
      
      # take a weighted mean of these distances. there should be no NA distances
      # except for the distance between a species and itself; set na.rm to TRUE
      indexGenResult[j] <- weighted.mean(theGenDist, theWeight, na.rm=TRUE)
      indexMorphResult[j] <- weighted.mean(theMorphDist, theWeight, na.rm=TRUE)
      
      # grab the abundance of spp[i] in that index
      tempAbund <- otherSpp$n[otherSpp$name==spp[i]]
      
      # calculate a weight as the proportion of lists in that index it occurred in
      # (its "abundance") divided by the max number of lists any species in the index
      # occurred in
      spWeight[j] <- tempAbund/max(otherSpp$n)
    }
  }
  
  # calculate the weighted mean distance to all other species in each of the indices
  # a species occurs in. use na.rm because some of the indices can be NA
  spGenResult <- weighted.mean(indexGenResult, spWeight, na.rm=TRUE)
  spMorphResult <- weighted.mean(indexMorphResult, spWeight, na.rm=TRUE)
  
  oneResult <- data.frame(spGenResult,spMorphResult)
  oneResult
})

finalResults <- Reduce(rbind, spResults)
finalResults <- data.frame(species=spp, finalResults)

# looks great. save out
write.csv(finalResults, "data/ebirdPhyloMorpho_3Mar2025.csv", row.names=FALSE)

