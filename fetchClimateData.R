rm(list=ls(all=TRUE)); graphics.off();
setwd("C:/cloud/MEGA/Projects/sApropos/")
require(RFc)
require(dplyr)
require(magrittr)
# Precipitation rate variable in FC is prate, units mm/month
# Temperature variable is airt, units degrees C


# read data ----------------------------------------------------------------------
d           <- read.csv("all_demog_6tr.csv", stringsAsFactors = F)

# exclude spp with no Lat/Lon info, AND no start/end year
d           <- subset(d, !is.na(Lat) & !is.na(Lon) & 
                         !is.na(MatrixStartYear) & !is.na(MatrixEndYear) )

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
                           clim_var = as.numeric(fc_raw$values[1,]),
                           stringsAsFactors = F)
  
  # change name of variable 
  fc_out    <- rename_(fc_out, .dots = setNames("clim_var", var) )
  
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


# download data separately (computer crashes otherwise) --------------------------------------------------

# air temperature data
spp_airt  <- lapply(1:10, climate_spp, "airt", grouped_data, 49)
airt_1    <- Reduce(function(...) rbind(...), spp_airt)
write.csv(airt_1, "C:/cloud/MEGA/Projects/sApropos/airt_fc_1.csv", row.names = F)
rm(airt_1)

spp_airt  <- lapply(11:18, climate_spp, "airt", grouped_data, 49)
airt_2    <- Reduce(function(...) rbind(...), spp_airt)
write.csv(airt_2, "C:/cloud/MEGA/Projects/sApropos/airt_fc_2.csv", row.names = F)
rm(airt_2)

spp_airt   <- lapply(20:30, climate_spp, "airt", grouped_data, 49)
airt_3    <- Reduce(function(...) rbind(...), spp_airt)
write.csv(airt_3, "C:/cloud/MEGA/Projects/sApropos/airt_fc_3.csv", row.names = F)
rm(airt_3)

spp_airt  <- lapply(31:40, climate_spp, "airt", grouped_data, 49)
airt_4    <- Reduce(function(...) rbind(...), spp_airt)
write.csv(airt_4, "C:/cloud/MEGA/Projects/sApropos/airt_fc_4.csv", row.names = F)
rm(airt_4)

spp_airt  <- lapply(41:50, climate_spp, "airt", grouped_data, 49)
airt_5    <- Reduce(function(...) rbind(...), spp_airt)
write.csv(airt_5, "C:/cloud/MEGA/Projects/sApropos/airt_fc_5.csv", row.names = F)
rm(airt_5)

spp_airt  <- lapply(51:60, climate_spp, "airt", grouped_data, 49)
airt_6    <- Reduce(function(...) rbind(...), spp_airt)
write.csv(airt_6, "C:/cloud/MEGA/Projects/sApropos/airt_fc_6.csv", row.names = F)
rm(airt_6)

spp_airt  <- lapply(61:70, climate_spp, "airt", grouped_data, 49)
airt_7    <- Reduce(function(...) rbind(...), spp_airt)
write.csv(airt_7, "C:/cloud/MEGA/Projects/sApropos/airt_fc_7.csv", row.names = F)
rm(airt_7)

spp_airt  <- lapply(71:78, climate_spp, "airt", grouped_data, 49)
airt_8    <- Reduce(function(...) rbind(...), spp_airt)
write.csv(airt_8, "C:/cloud/MEGA/Projects/sApropos/airt_fc_8.csv", row.names = F)
rm(airt_8)


# potential evapotranspiration data 
spp_prec  <- lapply(1:10, climate_spp, "prate", grouped_data, 49)
precip_1  <- Reduce(function(...) rbind(...), spp_prec)
write.csv(precip_1, "C:/cloud/MEGA/Projects/sApropos/precip_fc_1.csv", row.names = F)
rm(precip_1)

spp_prec  <- lapply(11:20, climate_spp, "prate", grouped_data, 49)
precip_2  <- Reduce(function(...) rbind(...), spp_prec)
write.csv(precip_2, "C:/cloud/MEGA/Projects/sApropos/precip_fc_2.csv", row.names = F)
rm(precip_2)

spp_prec   <- lapply(21:30, climate_spp, "prate", grouped_data, 49)
precip_3     <- Reduce(function(...) rbind(...), spp_prec)
write.csv(precip_3, "C:/cloud/MEGA/Projects/sApropos/precip_fc_3.csv", row.names = F)
rm(precip_3)

spp_prec   <- lapply(31:40, climate_spp, "prate", grouped_data, 49)
precip_4     <- Reduce(function(...) rbind(...), spp_prec)
write.csv(precip_4, "C:/cloud/MEGA/Projects/sApropos/precip_fc_4.csv", row.names = F)
rm(precip_4)

spp_prec   <- lapply(41:50, climate_spp, "prate", grouped_data, 49)
precip_5     <- Reduce(function(...) rbind(...), spp_prec)
write.csv(precip_5, "C:/cloud/MEGA/Projects/sApropos/precip_fc_5.csv", row.names = F)
rm(precip_5)

spp_prec   <- lapply(51:60, climate_spp, "prate", grouped_data, 49)
precip_6     <- Reduce(function(...) rbind(...), spp_prec)
write.csv(precip_6, "C:/cloud/MEGA/Projects/sApropos/precip_fc_6.csv", row.names = F)
rm(precip_6)

spp_prec   <- lapply(61:70, climate_spp, "prate", grouped_data, 49)
precip_7     <- Reduce(function(...) rbind(...), spp_prec)
write.csv(precip_7, "C:/cloud/MEGA/Projects/sApropos/precip_fc_7.csv", row.names = F)
rm(precip_7)

spp_prec   <- lapply(71:78, climate_spp, "prate", grouped_data, 49)
precip_8     <- Reduce(function(...) rbind(...), spp_prec)
write.csv(precip_8, "C:/cloud/MEGA/Projects/sApropos/precip_fc_8.csv", row.names = F)
rm(precip_8)


# potential evapotranspiration data 
spp_pet   <- lapply(1:10, climate_spp, "pet", grouped_data, 49)
pet_1     <- Reduce(function(...) rbind(...), spp_pet)
write.csv(pet_1, "C:/cloud/MEGA/Projects/sApropos/pet_fc_1.csv", row.names = F)
rm(pet_1)

spp_pet   <- lapply(11:20, climate_spp, "pet", grouped_data, 49)
pet_2     <- Reduce(function(...) rbind(...), spp_pet)
write.csv(pet_2, "C:/cloud/MEGA/Projects/sApropos/pet_fc_2.csv", row.names = F)
rm(pet_2)

spp_pet   <- lapply(21:30, climate_spp, "pet", grouped_data, 49)
pet_3     <- Reduce(function(...) rbind(...), spp_pet)
write.csv(pet_3, "C:/cloud/MEGA/Projects/sApropos/pet_fc_3.csv", row.names = F)
rm(pet_3)

spp_pet   <- lapply(31:40, climate_spp, "pet", grouped_data, 49)
pet_4     <- Reduce(function(...) rbind(...), spp_pet)
write.csv(pet_4, "C:/cloud/MEGA/Projects/sApropos/pet_fc_4.csv", row.names = F)
rm(pet_4)

spp_pet   <- lapply(41:50, climate_spp, "pet", grouped_data, 49)
pet_5     <- Reduce(function(...) rbind(...), spp_pet)
write.csv(pet_5, "C:/cloud/MEGA/Projects/sApropos/pet_fc_5.csv", row.names = F)
rm(pet_5)

spp_pet   <- lapply(51:60, climate_spp, "pet", grouped_data, 49)
pet_6     <- Reduce(function(...) rbind(...), spp_pet)
write.csv(pet_6, "C:/cloud/MEGA/Projects/sApropos/pet_fc_6.csv", row.names = F)
rm(pet_6)

spp_pet   <- lapply(61:70, climate_spp, "pet", grouped_data, 49)
pet_7     <- Reduce(function(...) rbind(...), spp_pet)
write.csv(pet_7, "C:/cloud/MEGA/Projects/sApropos/pet_fc_7.csv", row.names = F)
rm(pet_7)

spp_pet   <- lapply(71:78, climate_spp, "pet", grouped_data, 49)
pet_8     <- Reduce(function(...) rbind(...), spp_pet)
write.csv(pet_8, "C:/cloud/MEGA/Projects/sApropos/pet_fc_8.csv", row.names = F)
rm(pet_8)

