# Formatting functions for climate and lambda data-------------------------------------

# format species 
format_species <- function(spp_name, lam){
  
  lam_sel <- lam %>%
                subset(SpeciesAuthor == spp_name) %>%
                dplyr::select(MatrixEndYear, MatrixEndMonth, MatrixPopulation, lambda) %>%
                setNames(c("year","month","population","lambda")) %>%
                mutate(log_lambda = log(lambda),
                       population = as.factor(population) )
  
  # list of pop-specific lambdas
  lam_l   <- lam_sel %>%
                dplyr::select(-population) %>%
                split( lam_sel$population )
  
  return(lam_l)
  
}

# separate climate variables by population 
clim_list <- function(spp_name, clim, lam_spp){
  
  # select species-specific climate
  clim_spp  <- clim %>% 
                  subset(species == spp_name) %>%
                  mutate( population = as.factor(population) )
    
  # list of climate
  clim_l    <- clim_spp %>%
                  dplyr::select(-population) %>%
                  split( clim_spp$population )

  return(clim_l)
  
}

# detrend population-level climate; put it in "long" form
clim_detrend <- function(clim_x){ #, pops

  # format day one
  day_one   <- as.Date( paste0("1/1/", first(clim_x$year) ), 
                        format="%d/%m/%Y") 

  # climate data
  clim_d    <- as.Date(1:nrow(clim_x), day_one-1) %>%
                as.character %>%
                as.data.frame(stringsAsFactors=F) %>%
                separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
                bind_cols(clim_x) %>%
                dplyr::select(-year,-day) %>%
                setNames( c("year", "month", "day", "species", "ppt") )
  
  # monthly climates
  clim_m   <- clim_d %>%
                group_by(year, month) %>%
                summarise( ppt = sum(ppt, na.rm=T) )  %>%
                spread(month, ppt) %>%
                setNames( c("year",month.abb)) %>%            
                as.data.frame() #%>%
                #add_column(population = pops, .after = 1)
           
  # detrend climate
  d_prec   <- apply(clim_m[,-1], 2, FUN = scale, center = T, scale = T) %>%
                as.data.frame() %>%
                bind_cols(clim_m[,"year",drop=F]) %>%
                dplyr::select( c("year", month.abb) )

  return(d_prec)
  
}

# climate in long form 
clim_long <- function(clim_detr, lambda_data, m_back){
  
  # fecth year range, observation month
  years     <- lambda_data$year %>% unique
  yr_range  <- range(lambda_data$year)
  m_obs     <- lambda_data$month %>% unique
  
  # detrended climate in "long" form
  long_out  <- clim_detr %>%
                  subset(year < (yr_range[2]+1) & year > (yr_range[1] - 4) ) %>%
                  gather(month, precip, Jan:Dec) %>%
                  setNames(c("year", "month", "clim_value")) %>% 
                  mutate(month_num = factor(month, levels = month.abb) ) %>% 
                  mutate(month_num = as.numeric(month_num)) %>% 
                  arrange(year, month_num)
  
  # select temporal extent
  clim_back <- function(yrz, m_obs, dat) {
    id <- which(dat$year == yrz & dat$month_num == m_obs)
    r  <- c( id:(id - (m_back-1)) )
    return(dat[r,"clim_value"])
  }

  # climate data in matrix form 
  mat_form<- function(dat, years){
    do.call(rbind, dat) %>% 
      as.data.frame %>%
      tibble::add_column(year = years, .before=1)
  }
  
  # calculate monthly precipitation values
  precip_l  <- lapply(years, clim_back, m_obs, long_out)
  x_precip  <- mat_form(precip_l, years)
  return(x_precip)
  
}

# combine climate data frames (if any)
lambda_plus_clim <- function(lambdas_l, clim_mat_l){
  
  # lambda and climate n. of populations correspond?
  if( length(lambdas_l) != length(clim_mat_l) ) stop("lambda and climate lists have differing lengths")
  
  # add population name to data frames
  population_add <- function(x, pop_name){
    tibble::add_column(x, population = pop_name)
  }
  lambdas_l   <- Map(population_add, lambdas_l, names(lambdas_l) )
  clim_mat_l  <- Map(population_add, clim_mat_l, names(clim_mat_l) )
  
  # merge 
  if( length(lambdas_l) > 1){ # if n. of populations exceeds 1
    
    lambdas   <- Reduce(function(...) rbind(...), lambdas_l)
    climates  <- Reduce(function(...) rbind(...), clim_mat_l)
    clim_lam  <- merge(lambdas, climates)
    
  } else {
    
    clim_lam  <- merge(lambdas_l[[1]], clim_mat_l[[1]]) 
    
  }
    
  clim_lam    <- arrange(clim_lam, year, population)
  lam_out     <- dplyr::select(clim_lam, year:log_lambda)
  clim_out    <- dplyr::select(clim_lam, -c(year:log_lambda) )
  out         <- list(lambdas = lam_out, climate = clim_out)
  
  return(out)
  
}


# test functions ------------------------------------------------
# 
# setwd("C:/cloud/Dropbox/sAPROPOS project/DemogData")
# library(dplyr)
# library(tidyr)

# # read data -----------------------------------------------------
# lam     <- read.csv("lambdas_6tr.csv", stringsAsFactors = F) %>%
#             subset( !is.na(MatrixEndMonth) )
# clim    <- read.csv("precip_fc_demo.csv",  stringsAsFactors = F)
# spp     <- clim$species %>% unique


# read data -----------------------------------------------------

# spp_name      <- spp[3] # test run w/ spp number 1
# m_back        <- 24     # months back
# 
# # lambda data
# spp_lambdas   <- format_species(spp_name, lam)
# 
# # climate data
# clim_separate <- clim_list(spp_name, clim, spp_lambdas)
# clim_detrnded <- lapply(clim_separate, clim_detrend)
# clim_mats     <- Map(clim_long, clim_detrnded, spp_lambdas, m_back)
# 
# # model data
# mod_data  <- lambda_plus_clim(spp_lambdas, clim_mats)
