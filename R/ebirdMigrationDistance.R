library(data.table)
library(foreach)
library(doParallel)
library(h3r)
library(dplyr)
library(purrr)
library(clootl)
library(phytools)
library(viridis)

#set the working directory. this should be the directory that contains the data
#and R folders. if you got here by double-clicking the R script, the following
#will probably work
setwd("../")

# load all the ebird data

# RAM and compute issues very low so far, let's try to parallelize over these.
# read in all of your eBird parsed data. this is in 2023 taxonomy
#identify all the files/species you will loop over
allFiles <- list.files("data/queried&processed/")
spp <- unlist(lapply(strsplit(allFiles, "_"), "[", 1))

#append the correct file prefix for loading quickly
allFiles <- paste("data/queried&processed/", allFiles, sep="")

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
# combinations, into just index-month takes 54s
system.time(concat <- foreach(i = 1:length(results)) %dopar%
              {
                if(is.null(dim(results[[i]])))
                {
                  toStore <- NA
                }
                else
                {
                  # create the new index
                  tempIndex <- lapply(strsplit(results[[i]]$index, "-"), "[", 1)
                  tempMonth <- lapply(strsplit(results[[i]]$index, "-"), "[", 3)
                  results[[i]]$combined.index <- paste(tempIndex, tempMonth, sep="-")
                  results[[i]]$just.index <- unlist(tempIndex)
                  results[[i]]$month <- unlist(tempMonth)
                  
                  # group and summarize (summing abundances)
                  grouped <- group_by(results[[i]], combined.index)
                  toStore <- as.data.frame(summarize(grouped, n=sum(n),
                                                     month=unique(month),
                                                     just.index=unique(just.index)))
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
combinedDF <- imap_dfr(concat, ~ mutate(.x, species.name = .y))

# remove concat too
rm(concat)

# summarize the total number of individuals per index
grouped <- group_by(combinedDF, combined.index)
totalIndiv <- as.data.frame(summarize(grouped, total.indiv=sum(n)))

# merge in
combinedDF <- merge(combinedDF, totalIndiv, by="combined.index")

# remove some things for memory
rm(grouped)

# calculate a frequency from this
combinedDF$freq <- combinedDF$n/combinedDF$total.indiv

# figure out which species you will calculate this for. include only those species
# which have been observed in 10 months of the year
grouped <- group_by(combinedDF, species.name)
totalMonths <- as.data.frame(summarize(grouped, months=length(unique(month))))
spp <- totalMonths$species.name[totalMonths$months >= 10]

# now cut to these species
combinedDF <- combinedDF[combinedDF$species.name %in% spp,]

# find all unique grid cells
cells <- unique(combinedDF$just.index)

# find the lat and long of these
geolocs <- cellToLatLng(cells)
geolocs <- data.frame(index=cells, geolocs)

# merge these geo locations back in
combinedDF <- merge(combinedDF, geolocs, by.x="just.index", by.y="index")

# again remove some things
rm(grouped)

# loop over species and calculate a frequency-weighted latitude and longitude
# centroid per month
centroids <- list()

# find the frequency-weighted average geographic centroid using a function you
# worked up with ChatGPT!
average_geographic_centroid <- function(lat, lng, weights = NULL) {
  # Convert degrees to radians
  lat <- lat * pi / 180
  lng <- lng * pi / 180
  
  # Convert to Cartesian coordinates
  x <- cos(lat) * cos(lng)
  y <- cos(lat) * sin(lng)
  z <- sin(lat)
  
  # Compute weighted averages
  if (is.null(weights)) {
    # If no weights provided, use unweighted mean
    x_avg <- mean(x)
    y_avg <- mean(y)
    z_avg <- mean(z)
  } else {
    # Use weighted mean
    x_avg <- weighted.mean(x, weights)
    y_avg <- weighted.mean(y, weights)
    z_avg <- weighted.mean(z, weights)
  }
  
  # Convert back to latitude and longitude
  avg_lng <- atan2(y_avg, x_avg) * 180 / pi  # Convert to degrees
  avg_lat <- asin(z_avg / sqrt(x_avg^2 + y_avg^2 + z_avg^2)) * 180 / pi  # Convert to degrees
  
  return(c(avg_lat, avg_lng))
}

# now loop over species and actually implement
for(i in 1:length(spp))
{
  print(i)
  
  # subset to obs for that species
  temp <- combinedDF[combinedDF$species.name==spp[i],]
  
  # set aside the unique months
  months <- unique(temp$month)
  
  # set up a results frame
  output <- data.frame(month=months, latitude.centroid=0, longitude.centroid=0)
  
  # loop over the months now
  for(j in 1:length(months))
  {
    # subset temp to just the month
    temp2 <- temp[temp$month==months[j],]
    
    # calculate a weighted average using the new fxn
    tempResult <- average_geographic_centroid(lat=temp2$lat, lng=temp2$lng, weights=temp2$freq)
    output[j,"latitude.centroid"] <- tempResult[1]
    output[j,"longitude.centroid"] <- tempResult[2]
  }
  
  # set result into place
  centroids[[i]] <- output
}

# define a haversine function. note that this is a ChatGPT modified version of
# a haversine function you have used extensively.
haversine <- function(lon1, lat1, lon2, lat2) {
  
  if(!is.numeric(c(lon1, lat1, lon2, lat2)))
    stop("Inputs are not numeric")
  
  # Convert degrees to radians
  lon1 <- lon1 * pi / 180
  lat1 <- lat1 * pi / 180
  lon2 <- lon2 * pi / 180
  lat2 <- lat2 * pi / 180
  
  R <- 6371 # Earth mean radius [km]
  delta.lon <- (lon2 - lon1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) *
    sin(delta.lon/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  
  return(d) # Distance in km
}


# define a function that will sort a dataframe by month (tested that it sorts right
# even though month is stored as a charcter), duplicate the first row at the bottom,
# then calculate the haversine distance between every row and the next
sortAndStuff <- function(dat.frame)
{
  # sort the frame
  dat.frame <- dat.frame[order(dat.frame$month),]
  
  # find the last row
  lastRow <- dim(dat.frame)[1]
  
  # add the first as a new last
  dat.frame[1+lastRow,] <- dat.frame[1,]
  
  # set up a results vector and loop over frame. skip last row
  results <- c()
  for(i in 1:dim(dat.frame)[1]-1)
  {
    results[i] <- haversine(lat1=dat.frame[i,"latitude.centroid"],
                            lat2=dat.frame[i+1,"latitude.centroid"],
                            lon1=dat.frame[i,"longitude.centroid"],
                            lon2=dat.frame[i+1,"longitude.centroid"])
  }
  results
}

# use the function!
system.time(distances <- mclapply(centroids, sortAndStuff, mc.cores=6))

# summarize and give names
finalDists <- unlist(lapply(distances, sum))
names(finalDists) <- spp

# save out the results
writeOut <- data.frame(species.code=names(finalDists), dist=finalDists)
#write.csv(writeOut, "data/eBirdMigDist_1Dec2024.csv", row.names=FALSE)

