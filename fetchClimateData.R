rm(list=ls(all=TRUE)); graphics.off();
setwd("~/Dropbox/sAPROPOS project/DemogData")
require(RFc)
require(dplyr)
require(magrittr)
# Precipitation rate variable in FC is prate, units mm/month
# Temperature variable is airt, units degrees C


# read data --------------------------------------------------------
d           <- read.csv("lambdas_6tr.csv", stringsAsFactors = F)

# exclude spp with no Lat/Lon info, AND no start/end year
d           <- subset(d, !is.na(Lat) & !is.na(Lon) & 
                         !is.na(MatrixStartYear) & !is.na(MatrixEndYear))


# group data based on Species and population replicate ---------------------------
grouped_data  <- d %>%
                  group_by(SpeciesAuthor, MatrixPopulation) %>%
                  summarise(lat = first(Lat),
                            lon = first(Lon),
                            end_year = max(MatrixEndYear) )
                

# (nested) functions to fecth climate --------------------------------------------

# fetch daily climate
fetch_daily_clim <- function(yr, var, sp){

  fc_raw    <- fcTimeSeriesDaily(variable = var,
                                  latitude = sp$lat, longitude = sp$lon,
                                  firstYear = yr, lastYear = yr)
        
  # format data into a data frame
  fc_out    <- data.frame(species = as.character(sp$SpeciesAuthor),
                           population = sp$MatrixPopulation,
                           year = as.integer(yr), 
                           day = as.integer(fc_raw$days), 
                           ppt = as.numeric(fc_raw$values[1,]),
                           stringsAsFactors = F)
  
  return(fc_out)
  
}

# fetch climate across species
climate_spp <- function(sp_i, var, grouped_data, yr_back){
  
  sp      <- grouped_data[sp_i,]
  yr_r    <- seq(sp$end_year-yr_back, sp$end_year, by = 1)
  tmp     <- lapply(yr_r, fetch_daily_clim, var, sp)
  sp_clim <- Reduce(function(...) rbind(...), tmp) %>%
              as.data.frame(stringsAsFactors = F)
  return(sp_clim)
  
}


# download data ----------------------------------------------------------

# precipitation data
spp_prate <- lapply(1:nrow(grouped_data), climate_spp, "prate", grouped_data, 49)
precip    <- Reduce(function(...) rbind(...), spp_prate)

# temperature data
spp_t     <- lapply(1:nrow(grouped_data), climate_spp, "airt", grouped_data)
temp      <- Reduce(function(...) rbind(...), spp_t)


# store -------------------------------------------------------------------
write.csv(precip, "precip_fc.csv", row.names = F)
write.csv(temp,   "temp_fc.csv", row.names = F)
