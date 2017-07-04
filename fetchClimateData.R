rm(list=ls(all=TRUE)); graphics.off();
setwd("C:/Users/Aldo/Dropbox/sAPROPOS project/DemogData")
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

# separate species: one spatial replicate/multiple spatial replicates
spat_rep    <- d %>%
                  select(SpeciesAuthor, Lat, Lon) %>%
                  unique %>%
                  group_by(SpeciesAuthor) %>%
                  summarise(spatial_rep = n() )
unrep_spp   <- subset(spat_rep, spatial_rep == 1)$SpeciesAuthor
rep_spp     <- subset(spat_rep, spatial_rep != 1)$SpeciesAuthor


# one spatial replicate only ---------------------------------------
unrep_d       <- subset(d, SpeciesAuthor %in% unrep_spp)

# Get per-study lat/lon and year range
grouped_data  <- unrep_d %>% 
                  group_by(SpeciesAuthor) %>%
                  summarise(lat = first(Lat), 
                            lon = first(Lon),
                            end_year = max(MatrixEndYear) ) 
                  

# (nested) functions to fecth climate --------------------------------------------

# fetch daily climate
fetch_daily_clim <- function(yr, var, sp){

  tmp   <- fcTimeSeriesDaily(variable = var,
                            latitude = sp$lat, longitude = sp$lon,
                            firstYear = yr, lastYear = yr)
  
  # format data into a data frame
  df_out <- data.frame(species = as.character(sp$SpeciesAuthor),
                       year = as.integer(yr), 
                       day = as.integer(tmp$days), 
                       ppt = as.numeric(tmp$values[1,]),
                       stringsAsFactors = F)
  
  return(df_out)
  
}

# fetch climate across species
climate_spp <- function(sp_i, var, grouped_datam, yr_back){
  
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


# store -------------------------------------------------------------
write.csv(precip, "precip_fc.csv", row.names = F)
write.csv(temp,   "temp_fc.csv", row.names = F)
