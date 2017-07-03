rm(list=ls(all=TRUE)); graphics.off();
# Uses FetchClimate package from CRAN
require(RFc)
require(dplyr)


# Precipitation rate variable in FC is prate, units mm/month
# Temperature variable is airt, units degrees C
setwd("C:/Users/Aldo/Dropbox/sAPROPOS project/DemogData")

# read data
data <- read.csv("vitalRates.csv", head=T)

# more than one location per species?
# d <- read.csv("lambdas_6tr.csv", head=T)
# grouped_data <- d %>% 
#                   select(SpeciesAuthor, Lat, Lon) %>%
#                   unique %>%
#                   group_by(SpeciesAuthor) %>%
#                   summarise(n = n() )


# Get per-study lat/lon and year range
grouped_data <- data %>% 
                  group_by(SpeciesAuthor) %>%
                  summarise(lat = first(Lat), 
                            lon = first(Lon),
                            firstYear = min(MatrixStartYear), 
                            lastYear = max(MatrixEndYear) )


# Use daily time series to get each year precipitation data ----------------------
ppt_daily_data <- data.frame(species = character(),
                             year=integer(),
                             day=integer(),
                             ppt=double(),
                             stringsAsFactors=FALSE)

for (i in seq_len(nrow(grouped_data))) {
  sp <- grouped_data[i,]
  for (yr in seq(from = sp$firstYear-2, to = sp$lastYear)) {
    tmp <- fcTimeSeriesDaily(variable="prate",
                             latitude = sp$lat, longitude = sp$lon,
                             firstYear = yr, lastYear = yr)
    tmp <- cbind(species = as.character(sp$SpeciesAuthor),
                 year = yr, day = tmp$days, ppt = tmp$values[1,])
    ppt_daily_data <- rbind(ppt_daily_data, tmp)
  }
}

write.csv(ppt_daily_data, "fcPpt.csv")


# Use daily time series to get each year temp data -----------------------------
temp_daily_data <- data.frame(species = character(),
                              year=integer(),
                              day=integer(),
                              temp=double(),
                              stringsAsFactors=FALSE)
for (i in seq_len(nrow(grouped_data))) {
  sp <- grouped_data[i,]

  for (yr in seq(from = sp$firstYear-2, to = sp$lastYear)) {
    tmp <- fcTimeSeriesDaily(variable="airt",
                             latitude = sp$lat, longitude = sp$lon,
                             firstYear = yr, lastYear = yr)
    tmp <- cbind(species = as.character(sp$SpeciesAuthor),
                 year = yr, day = tmp$days, temp = tmp$values[1,])
    temp_daily_data <- rbind(temp_daily_data, tmp)
  }
}

write.csv(temp_daily_data, paste0(wd, "/DemogData/fcTemp.csv"))
