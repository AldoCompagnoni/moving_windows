# species-specific lambda and climate data 
setwd("C:/cloud/Dropbox/sAPROPOS project/DemogData")
library(dplyr)
library(tidyr)


# read data -----------------------------------------------------
lam     <- read.csv("lambdas_6tr.csv", stringsAsFactors = F)
clim    <- read.csv("precip_fc.csv",  stringsAsFactors = F) #clim<-read.csv("precip_fc_demo.csv",  stringsAsFactors = F)
spp     <- clim$species %>% unique
m_back  <- 24


# format species ------------------------------------------------
format_species <- function(spp_name, lam){
  
  lam   <- lam %>%
            subset(SpeciesAuthor == spp_name) %>%
            dplyr::select(MatrixEndYear, MatrixPopulation, lambda) %>%
            setNames(c("year","MatrixPopulation","lambda")) %>%
            mutate(log_lambda = log(lambda))
  return(lam)
  
}

# format climate 
format_climate <- function(spp_name, clim, lam_spp){
  
  # select species-specific climate
  clim_spp  <- clim %>% 
                  subset(species == spp_name) %>%
                  mutate( population = as.factor(population) )
    
  # list of climate
  clim_l    <- clim_spp %>%
                  select(-population) %>%
                  split( clim_spp$population )

  return(clim_l)
  
}

# format populations
format_pop <- function(clim_x){

  # format day one
  day_one   <- as.Date( paste0("1/1/", first(clim_x$year) ), 
                        format="%d/%m/%Y") 

  # climate data
  clim_d    <- as.Date(1:nrow(clim_x), day_one-1) %>%
                as.character %>%
                as.data.frame(stringsAsFactors=F) %>%
                separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
                bind_cols(clim_x) %>%
                select(-year,-day) %>%
                setNames( c("year", "month", "day", "species", "ppt") )
  
  # monthly climates
  clim_m   <- clim_d %>%
                group_by(year, month) %>%
                summarise( ppt = sum(ppt, na.rm=T) )  %>%
                spread(month, ppt) %>%
                setNames( c("year",month.abb))
                
  # detrend climate
  d_prec   <- apply(clim_m[,-1], 2, FUN = scale, center = T, scale = T) %>% 
                as.data.frame() %>%
                bind_cols(clim_m[,1]) %>%
                select( c("year",month.abb) )

  return(d_prec)
  
}


# format by species ----------------------------------------------------------
spp_name      <- spp[1] # test run w/ spp number 1

spp_lambdas   <- format_species(spp_name, lam)
spp_raw_clim  <- format_climate(spp_name, clim, lam_spp)
spp_clim      <- lapply(spp_raw_clim, format_pop)
