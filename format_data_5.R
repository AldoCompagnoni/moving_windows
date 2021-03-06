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
                  subset( species == spp_name ) %>%
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
clim_detrend <- function(clim_x, clim_var = "precip", foo){ 
  
  # precipitation can only be calculated via sums
  if( clim_var == "precip") foo = sum
  
  # years
  years_vec <- clim_x$year %>% unique %>% sort
  
  # clim_x value rename
  clim_x    <- clim_x %>% setNames( c("species","year","day","value") )
  
  # prepare indexes
  split_id  <- Filter(function(x) x < 366, which( (0:365 %% 5) == 0 ) ) %>%
                  lapply( function(x) x + (0:4) ) %>% 
                  setNames( paste0(c(1:length(.))) )
  
  # period summaries
  period_summ <- function(ii,clim_yr, foo){
    
    #fun         <- substitute( foo )
    period_subs <- clim_yr[ii,]
    period_clim <- data.frame( species = unique(period_subs$species), 
                               year    = unique(period_subs$year), 
                               value   = foo(period_subs$value, na.rm=T) %>% eval
                              )
    return(period_clim)
    
  }
  
  # repeat summaries by year
  years_pop <- function(yrs, clim_x, period_summ, split_id, foo){
    
    clim_yr    <- subset(clim_x, year == yrs)
    # create summaries data frame
    period_l  <- lapply(split_id, period_summ, clim_yr, foo)
    period_df <- Reduce(function(...) rbind(...), period_l) %>% 
                    tibble::add_column( period = 1:length(split_id), .after=2 ) %>%
                    dplyr::select( -species ) %>%
                    spread( period, value )
    return(period_df)
    
  }
  
  # substitute -Inf with NA
  na_for_inf  <- function(x) replace(x, x==-Inf, NA)
  
  # 5 day interval summaries
  clim_int_l  <- lapply(years_vec, years_pop, clim_x, period_summ, split_id, foo)
  clim_int    <- Reduce(function(...) rbind(...), clim_int_l) %>%
                    apply(2, na_for_inf)
  
  # detrend climate 
  d_prec      <- apply(clim_int[,-1], 2, FUN = scale, center = T, scale = T) %>%
                    as.data.frame() %>%
                    bind_cols( as.data.frame(clim_int[,"year", drop=F]) ) %>%
                    dplyr::select( c("year", 1:(ncol(clim_int)-1) ) )
  
  # Present 
  nans_n <- d_prec[,-1] %>% as.matrix %>% is.na %>% sum
  if(nans_n > 0) print( "Warning: NANs present in climate predictor" )
  
  # Make NAs 0
  for(c_i in 1:ncol(d_prec) ){
    d_prec[,c_i] <- replace(d_prec[,c_i], is.na(d_prec[,c_i]), 0)
  }
  
  return(d_prec)
  
}

# climate in long form 
clim_long <- function(clim_detr, lambda_data, m_back){
  
  # fecth year range, observation month
  years     <- lambda_data$year %>% unique
  yr_range  <- range(lambda_data$year)
  m_obs     <- lambda_data$month %>% unique
  
  # format day one (just a trick)
  day_one   <- as.Date( paste0("1/1/", first(clim_detr$year) ), 
                        format="%d/%m/%Y" ) 
  
  # dates data frame (with julian date)
  dates_df  <- data.frame( date = as.Date(0:364, day_one) ) %>%
                  mutate( date  = as.character(date) ) %>%
                  mutate( jul   = 1:365 ) %>%
                  separate( date, into = c("year","month","day"), sep = "-") %>%
                  mutate( month = as.numeric(month) )
  
  # identify what's the 5-day period where month observation is located 
  five_day_id <- function(ii, x){
    # pick last day of the month
    out <- subset(x, month == ii) %>% .$jul %>% last  
    # number of 5 day period
    return( round(out / 5) ) 
  }
  
  # 5-day period of month (demographic) observation
  five_d_id   <- five_day_id(m_obs, dates_df)
  
  # N. of 5-day periods "back"
  five_d_back <- round( c(m_back * 30.41667) / 5)
  
  # detrended climate in "long" form
  long_out  <- clim_detr %>%
                  subset(year < (yr_range[2]+1) & year > (yr_range[1] - 6) ) %>%
                  gather(five_d, clim_value, paste0(1:73)) %>%
                  mutate(five_d = as.numeric(five_d) ) %>% 
                  arrange(year, five_d)
  
  # select temporal extent
  clim_back <- function(yrz, long_dat, five_d_id, five_d_back){
    id <- which(long_dat$year == yrz & long_dat$five_d == five_d_id)
    r  <- c( id:(id - (five_d_back-1)) )
    return(long_dat[r,"clim_value"])
  }
  clim_x_l  <- lapply(years, clim_back, long_out, five_d_id, five_d_back)
  
  # climate data in matrix form 
  year_by_month_mat <- function(dat, years){
    do.call(rbind, dat) %>% 
      as.data.frame %>%
      tibble::add_column(year = years, .before=1)
  }
  
  # calculate monthly precipitation values
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
  
  # variables
  var_select  <- c("year", "month", "population", response)
  
  # order, and erase cases in which lambda == 0 (e.g. Eryngium_alpinum, BOU, year 2009)
  clim_lam    <- arrange(clim_lam, year, population)
  # erase any row containing NAs (for Dalgleish et al. 2010 data)
  r_id        <- lapply(clim_lam, function(x) which(is.na(x)) ) %>% unlist
  if( length(r_id) > 0 ) clim_lam  <- clim_lam[-r_id,]
  
  resp_out    <- dplyr::select(clim_lam, var_select)
  clim_out    <- dplyr::select(clim_lam, -which( names(clim_lam) %in% var_select ) )
  out         <- list(resp = resp_out, climate = clim_out)
  
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
  