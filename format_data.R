# species-specific lambda and climate data 
setwd("C:/cloud/Dropbox/sAPROPOS project/DemogData")
library(dplyr)
library(tidyr)


# read data -----------------------------------------------------
lam     <- read.csv("lambdas_6tr.csv", stringsAsFactors = F)
clim    <- read.csv("precip_fc.csv",  stringsAsFactors = F)
spp     <- clim$species %>% unique
m_back  <- 24


# format species ------------------------------------------------
format_species <- function(spp_name, lam){
  
  lam   <- lam %>%
            subset(SpeciesAuthor == spp_name) %>%
            dplyr::select(MatrixEndYear, lambda) %>%
            setNames(c("year","lambda")) %>%
            mutate(log_lambda = log(lambda))
  return(lam)
  
}

lam_spp <- format_species(spp_name, lam)

# format climate ------------------------------------------------
format_climate <- function(spp_name, clim, lam_spp){
  
  # select species-specific climate
  clim_spp <- clim %>% subset(species == spp_name)
  
  # format day one
  day_one  <- as.Date(paste0("1/1/",first(clim_spp$year)), 
                      format="%d/%m/%Y") 
  
  # climate data
  clim_d   <- as.Date(1:nrow(clim_spp), day_one-1) %>%
                as.character %>%
                as.data.frame(stringsAsFactors=F) %>%
                separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
                bind_cols(clim_spp) %>%
                select(-year,-day) %>%
                setNames( c("year", "month", "day", "species", "ppt") )
       
  # monthly climates
  clim_m   <- clim_d %>%
                group_by(year, month) %>%
                summarise( ppt = sum(ppt) ) %>%
                spread(month, ppt) %>%
                setNames( c("year",month.abb))
   
  # detrend climate
  d_prec   <- apply(clim_m[,-1], 2, FUN = scale, center = T, scale = T) %>% 
                as.data.frame() %>%
                bind_cols(clim_m[,1]) %>%
                select( c("year",month.abb) )
    
  # long form (again!)
  
  
  
}



# select species
spp_dur <- lam %>% 
            group_by(SpeciesAccepted) %>% 
            summarise(duration = length(unique(MatrixStartYear)))

# species list
spp_list<- lam$SpeciesAuthor %>% unique
  subset(Lat == 38.8 & Lon == -99.2) %>%
  .[,"SpeciesAccepted"] %>%
  as.character %>%
  unique

