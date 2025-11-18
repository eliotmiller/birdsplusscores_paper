library(dplyr)
library(sf)
library(terra)
library(ebirdst)

#set the working directory. this should be the directory that contains the data
#and R folders. if you got here by double-clicking the R script, the following
#will probably work
setwd("../")

# set the key. use your own key
#set_ebirdst_access_key("")

# figure out what species are available
avail <- ebirdst_runs

# go ahead and download all species' data, the products you think you need,
# and at just 27km resolution

# curious to try new syntax I'm picking up on from Python.
# because force = FALSE, you can probably re-run this whenever without
# worrying about downloding it all again
spp <- avail$species_code

# for(species in spp)
# {
#   ebirdst_download_status(
#     species,
#     path = ebirdst_data_dir(),
#     download_abundance = TRUE,
#     download_occurrence = TRUE,
#     download_count = FALSE,
#     download_ranges = TRUE,
#     download_regional = FALSE,
#     download_pis = TRUE,
#     download_ppms = TRUE,
#     download_all = FALSE,
#     pattern = "_27km_",
#     dry_run = FALSE,
#     force = FALSE,
#     show_progress = TRUE
#   )
# }

# now process the large datasets into species-level measures
# using slightly modified code from matt strimas-mackey to calculate
# resident/breeding season abundance and range size
medianAbundance <- c()
rangeSize <- c()

for(i in 1:length(spp))
{
  print(spp[i])
  
  # figure out if the species is a resident or not
  if(avail$is_resident[avail$species_code==spp[i]])
  {
    abd_breeding <- load_raster(spp[i],
                                period = "seasonal",
                                resolution = "27km") |> 
      subset("resident") |> 
      # remove zeros
      subst(0, NA)
  }
  else
  {
    # at least one species only has a "non-breeding" category here. check for it
    # and handle accordingly
    abd_breeding <- load_raster(spp[i],
                                period = "seasonal",
                                resolution = "27km")
    
    if(sum(names(abd_breeding) %in% "breeding"))
    {
      abd_breeding <- subset(abd_breeding, "breeding") |> 
        # remove zeros
        subst(0, NA)
    }
    
    else
    {
      abd_breeding <- subset(abd_breeding, "nonbreeding") |> 
        # remove zeros
        subst(0, NA)
    }
  }
  
  temp <- as.data.frame(abd_breeding, cells = TRUE, xy = TRUE, na.rm = TRUE)
  rangeSize[i] <- dim(temp)[1]
  medianAbundance[i] <- as.numeric(global(abd_breeding, fun = median, na.rm = TRUE))
}

# create a results obj
statusResults <- data.frame(species=spp, median.abund=medianAbundance, range.size=rangeSize)

# download the trends data now. see what is available
trends_runs <- ebirdst_runs |> 
  filter(has_trends) |> 
  select(species_code, common_name,
         trends_season, trends_region,
         trends_start_year, trends_end_year,
         trends_start_date, trends_end_date,
         rsquared, beta0)
glimpse(trends_runs)

spp <- trends_runs$species_code

# download the actual data
for(species in spp)
{
  ebirdst_download_trends(species)
}

# calculate a mean abundance weighted trend per species
meanTrend <- c()

for(species in spp)
{
  # load the data
  spTrendDF <- load_trends(species)
  
  # calculate avearge abundance-weighted trend
  meanTrend[species] <- sum(spTrendDF$abd * spTrendDF$abd_ppy) / sum(spTrendDF$abd)
}

names(meanTrend) <- spp

# merge and save out. drop the example set
toMerge <- data.frame(species=spp, trend=meanTrend)
results <- merge(statusResults, toMerge, all.x=TRUE)
results <- results[results$species!="yebsap-example",]

write.csv(results, "data/processedST.csv", row.names=FALSE)

