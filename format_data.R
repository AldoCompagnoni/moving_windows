# Formatting functions for climate and lambda data-------------------------------------

# species from Dalgleish et al. 2010 (fetchClimate data not useful in this case)
Dalgleish_spp <- c("Cirsium_undulatum", "Echinacea_angustifolia", 
                   "Hedyotis_nigricans","Lesquerella_ovalifolia", 
                   "Paronychia_jamesii", "Psoralea_tenuiflora",      
                   "Ratibida_columnifera", "Solidago_mollis", 
                   "Sphaeralcea_coccinea", "Thelesperma_megapotamicum")

# format species 
format_species <- function(spp_name, lam, response = "lambda"){
  
  # only in the case of lambda, we include "log_lambda" as well
  # if( response == "lambda" ) response = c("lambda", "log_lambda")

  # fetch what you need from 'lam' object
  lam_sel <- lam %>%
                subset( SpeciesAuthor == spp_name ) %>%
                dplyr::select( c("MatrixEndYear", "MatrixEndMonth", "MatrixPopulation", response) ) %>%
                setNames( c("year","month","population", response) ) %>%
                mutate( population = as.factor(population) )

  # list of pop-specific lambdas
  lam_l   <- lam_sel %>%
                dplyr::select( -population ) %>%
                split( lam_sel$population )
  
  return(lam_l)
  
}

# separate climate variables by population 
clim_list <- function(spp_name, clim, lam_spp){ # clim_var, 
 
  clim_spp  <- clim %>%
                  subset( species == spp_name) %>%
                  mutate( population = as.factor(population) ) 
                  
  if( !(spp_name %in% Dalgleish_spp) ){
    clim_spp <- subset(clim_spp, year > 1947 )
  }
  
  # climate variables list
  clim_l    <- clim_spp %>%
                  dplyr::select(-population) %>%
                  split( clim_spp$population )
  
  return(clim_l)
  
}

# detrend population-level climate; put it in "long" form
clim_detrend <- function(clim_x, clim_var = "precip", st_dev = FALSE ){ #, pops

  # format day one
  day_one   <- as.Date( paste0("1/1/", first(clim_x$year) ), 
                        format="%d/%m/%Y" ) 

  # climate data
  clim_d    <- as.Date(1:nrow(clim_x), day_one-1) %>%
                  as.character %>%
                  as.data.frame(stringsAsFactors=F) %>%
                  separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
                  bind_cols(clim_x) %>%
                  dplyr::select(-year,-day) %>%
                  setNames( c("year", "month", "day", "species", "value") )
    
  # # if climate_var airt, then do means, otherwise, do sums! 
  if( clim_var == "airt"){
    clim_m  <- clim_d %>%
      group_by(year, month) %>%
      summarise( value = mean(value, na.rm=T) )  %>%
      spread( month, value ) %>%
      setNames( c("year",month.abb) ) %>%            
      as.data.frame()
  } else{
    clim_m   <- clim_d %>%
      group_by(year, month) %>%
      summarise( value = sum(value, na.rm=T) )  %>%
      spread( month, value ) %>%
      setNames( c("year",month.abb) ) %>%            
      as.data.frame()
  }

  # if st_dev == T, this overrides the above conditional statements
  if(st_dev == T){
    clim_m   <- clim_d %>%
      group_by(year, month) %>%
      summarise( value = sd(value, na.rm=T) )  %>%
      spread( month, value ) %>%
      setNames( c("year",month.abb) ) %>%            
      as.data.frame()
  }
  
  # throw error
  if( !any( clim_var %in% c("precip","pet","airt","gdd")) ) {
    stop( paste0(clim_var," is not a supported varible") ) 
  }
  
  # detrend climate - but NOT if you are using GDD 
  if( clim_var != "gdd" ){
    d_prec   <- apply(clim_m[,-1], 2, FUN = scale, center = T, scale = T) %>%
                  as.data.frame() %>%
                  bind_cols( clim_m[,"year", drop=F] ) %>%
                  dplyr::select( c("year", month.abb) )
  } else{
    d_prec   <- clim_m
  }
  
  # Present 
  nans_n <- d_prec[,-1] %>% as.matrix %>% is.nan %>% sum
  if(nans_n > 0) print( "Warning: NANs present in climate predictor" )
  
  # Make NaNs 0
  for(c_i in 1:ncol(d_prec) ){
    d_prec[,c_i] <- replace(d_prec[,c_i], is.nan(d_prec[,c_i]), 0)
  }
  
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
                  subset(year < (yr_range[2]+1) & year > (yr_range[1] - 6) ) %>%
                  gather(month, precip, Jan:Dec) %>%
                  setNames(c("year", "month", "clim_value") ) %>% 
                  mutate(month_num = factor(month, levels = month.abb) ) %>% 
                  mutate(month_num = as.numeric(month_num) ) %>% 
                  arrange(year, month_num)
  
  # select temporal extent
  clim_back <- function(yrz, m_obs, dat){
    id <- which(dat$year == yrz & dat$month_num == m_obs)
    r  <- c( id:(id - (m_back-1)) )
    return(dat[r,"clim_value"])
  }

  # climate data in matrix form 
  year_by_month_mat <- function(dat, years){
    do.call(rbind, dat) %>% 
      as.data.frame %>%
      tibble::add_column(year = years, .before=1)
  }
  
  # calculate monthly precipitation values
  clim_x_l  <- lapply(years, clim_back, m_obs, long_out)
  x_clim    <- year_by_month_mat(clim_x_l, years)
  return(x_clim)
  
}

# combine climate data frames (if any)
lambda_plus_clim <- function(lambdas_l, clim_mat_l, response = "lambda"){
  
  # lambda and climate n. of populations correspond?
  if( length(lambdas_l) != length(clim_mat_l) ) stop("response and climate lists have differing lengths")
  
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
    clim_lam  <- full_join(lambdas, climates)
    
  } else {
    
    clim_lam  <- full_join(lambdas_l[[1]], clim_mat_l[[1]]) 
    
  }
  
  # if( response == "lambda"){
  #   # order, and erase cases in which lambda == 0 (e.g. Eryngium_alpinum, BOU, year 2009)
  #   clim_lam    <- arrange(clim_lam, year, population)  %>%
  #                     subset( lambda != 0 ) 
  #   # erase any row containing NAs (for Dalgleish et al. 2010 data)
  #   r_id        <- lapply(clim_lam, function(x) which(is.na(x)) ) %>% unlist
  #   if( length(r_id) > 0 ) clim_lam  <- clim_lam[-r_id,]
  #   lam_out     <- dplyr::select(clim_lam, year:log_lambda)
  #   clim_out    <- dplyr::select(clim_lam, -c(year:log_lambda) )
  #   out         <- list(lambdas = lam_out, climate = clim_out)
  # }else{
  
  # order, and erase cases in which lambda == 0 (e.g. Eryngium_alpinum, BOU, year 2009)
  clim_lam    <- arrange(clim_lam, year, population)
  # erase any row containing NAs (for Dalgleish et al. 2010 data)
  r_id        <- lapply(clim_lam, function(x) which(is.na(x)) ) %>% unlist
  if( length(r_id) > 0 ) clim_lam  <- clim_lam[-r_id,]
  eval(parse(n=1, text=paste0("resp_out <- dplyr::select(clim_lam, year:",response,")")))
  eval(parse(n=1, text=paste0("clim_out <- dplyr::select(clim_lam, -c(year:",response,"))")))
  out         <- list(resp = resp_out, climate = clim_out)
  # }
  
  return(out)
  
}


# observed cliamtic range
observed_clim_range <- function(clim_x, lambda_d, spp_name, clim_var){

  # format day one
  day_one   <- as.Date( paste0("1/1/", first(clim_x$year) ), 
                        format="%d/%m/%Y" )
  
  # climate data
  clim_d    <- as.Date(1:nrow(clim_x), day_one-1) %>%
                  as.character %>%
                  as.data.frame(stringsAsFactors=F) %>%
                  separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
                  bind_cols(clim_x) %>%
                  dplyr::select(-year,-day) %>%
                  setNames( c("year", "month", "day", "species", "ppt") )
  
  # monthly climates
  # clim_var == "airt" calculate MEAN monthly air temperature 
  if( clim_var == "airt" ){
    clim_m   <- clim_d %>%
                    group_by(year, month) %>%
                    summarise( ppt = mean(ppt, na.rm=T) ) %>%
                    ungroup %>%
                    mutate( month = as.numeric(month) ) %>%
                    mutate( year = as.numeric(year) ) 
  }else{
    clim_m   <- clim_d %>%
                    group_by(year, month) %>%
                    summarise( ppt = sum(ppt, na.rm=T) ) %>%
                    ungroup %>%
                    mutate( month = as.numeric(month) ) %>%
                    mutate( year = as.numeric(year) ) 
  }
  
  # range of years
  max_yr    <- max(lambda_d$year)
  min_yr    <- min(lambda_d$year)
  month_i   <- unique(lambda_d$month)
  
  yearly_climate <- function(yrs){
    
    year_clim <- clim_m %>% subset( year == yrs & month < month_i + 1 )
    
    if( unique(lambda_d$month) != 12 ){
     
     append     <- clim_m %>% 
                      subset( year == (yrs-1) & month > month_i )
     year_clim  <- rbind(year_clim, append)
     
    }
  
    # Calculate means for "airt"
    if( clim_var == "airt"){
      out <- data.frame( year = yrs, ppt = mean(year_clim$ppt) )
    }else{
      out <- data.frame( year = yrs, ppt = sum(year_clim$ppt) )
    }
    
    return(out)
    
  }
  
  # yearly climates
  all_yrs   <- (max_yr-48):max_yr
  yr_clim_l <- lapply(all_yrs, yearly_climate)
  full_clim <- Reduce(function(...) rbind(...), yr_clim_l)
  obs_clim  <- subset(full_clim, year >= min_yr & year <= max_yr )
  
  # climate in full data set
  full_rng  <- max(full_clim$ppt) - min(full_clim$ppt)
  full_mean <- mean(full_clim$ppt)
  full_med  <- median(full_clim$ppt)
  full_sd   <- sd(full_clim$ppt)
  full_dev  <- range(abs(full_mean - full_clim$ppt))
  full_dev_r<- range(abs(full_med - full_clim$ppt))
  
  # observed range of climate anomalies
  obs_dev   <- range( abs(full_mean - obs_clim$ppt) )
  obs_dev_r <- range( abs(full_med - obs_clim$ppt) )
  obs_range <- max(obs_clim$ppt) - min(obs_clim$ppt)
  extr_yr_n <- sum(full_clim$ppt > max(obs_clim$ppt)) + sum(full_clim$ppt < min(obs_clim$ppt)) 
  obs_mean  <- mean(obs_clim$ppt)
  
  # proportions and means
  prop_yrs  <- (48-extr_yr_n) / 48 
  prop_rang <- obs_range / full_rng
  prop_var  <- (obs_dev[2] - obs_dev[1]) / (full_dev[2] - full_dev[1])
  prop_var_r<- (obs_dev_r[2] - obs_dev_r[1]) / (full_dev_r[2] - full_dev_r[1])
  mean_dev  <- (full_mean - obs_mean) / full_sd
  
  return( data.frame( species   = spp_name, 
                      prop_rang = prop_rang, 
                      prop_yrs  = prop_yrs,
                      prop_var  = prop_var,
                      prop_var_r= prop_var_r,
                      mean_dev  = mean_dev,
                      mean_clim = full_mean )
         )
  
}
  