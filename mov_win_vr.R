#bjtwd<-"C:/Users/admin_bjt162/Dropbox/A.Current/Ongoing_Collab_Research/sApropos project/"
rm(list=ls())
setwd("C:/cloud/Dropbox/sApropos/")
source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(dismo)
library(mgcv)
library(testthat)
library(rstan)
library(loo)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# "pipeable" Reduce rbind
rbind_l <- function(x) Reduce(function(...) rbind(...), x)

# climate predictor, response, months back, max. number of knots
response  <- "fec"
clim_var  <- "airt"
m_back    <- 36    
st_dev    <- FALSE

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv("all_demog_6tr.csv", stringsAsFactors = F)
m_info    <- read.csv("MatrixEndMonth_information.csv", stringsAsFactors = F)
clim      <- data.table::fread(paste0(clim_var,"_chelsa_hays.csv"),  stringsAsFactors = F)
# clim      <- data.table::fread(paste0(clim_var,"_fc_hays.csv"),  stringsAsFactors = F)

spp       <- lam$SpeciesAuthor %>% unique

# format data --------------------------------------------------------------------------------------

# set up model "family" based on response
if( response == "surv" | response == "grow" )             family = "beta" 
if( response == "fec" )                                   family = "gamma"
if( grepl("PreRep", response) | grepl("Rep", response) )  family = "beta"
if( response == "rho" | response == "react_fsa" )         family = "gamma"
if( response == "log_lambda")                             family = "normal"

expp_beta     <- 20

# set species (I pick Sphaeraclea_coccinea)
ii            <- 13
spp_name      <- spp[ii]

# lambda data
spp_resp      <- format_species(spp_name, lam, response)

# climate data
clim_separate <- clim_list(spp_name, clim, spp_resp)
#clim_detrnded <- lapply(clim_separate, clim_detrend, clim_var)
clim_detrnded <- lapply(clim_separate, clim_detrend, clim_var)
clim_mats     <- Map(clim_long, clim_detrnded, spp_resp, m_back)

# model data
mod_data          <- lambda_plus_clim(spp_resp, clim_mats, response)
mod_data$climate  <- mod_data$climate #/ diff(range(mod_data$climate))

# throw error if not enough data
if( nrow(mod_data$resp) < 6 ) stop( paste0("not enough temporal replication for '", 
                                              spp_name, "' and response variable '" , response, "'") )


# Transform response variables (if needed) ------------------------------------------------------------------

# transform survival/growth - ONLY if less than 30% data points are 1/0
if( response == "surv" | response == "grow" | grepl("PreRep", response) | grepl("Rep", response) ){
  
  raw_x <- mod_data$resp[,response]
  pc_1  <- sum( raw_x == 1 ) / length(raw_x)
  pc_0  <- sum( raw_x == 0 ) / length(raw_x)
  
  # for survival
  if( grepl("surv", response, ignore.case = T) & pc_1 < 0.3 ){
    n     <- length(raw_x)
    new_x <- ( raw_x*(n - 1) + 0.5 ) / n
    mod_data$resp[,response] <- new_x
  }
  
  # for growth
  if( grepl("grow", response, ignore.case = T) & pc_0 < 0.3 ){
    n     <- length(raw_x)
    new_x <- ( raw_x*(n - 1) + 0.5 ) / n
    mod_data$resp[,response] <- new_x
  }
  
}

# avoid absolute zeros
if( response == "fec" ){
  # transform from [0, infinity) to (0, infinity) 
  # I add quantity 2 orders of mag. lower than lowest obs value.
  mod_data$resp[,response] <- mod_data$resp[,response] + 1.54e-12 
} 

if( response == "rho" | response == "react_fsa" ){
  # bound responses to (0,infinity) instead of [1, infinity) 
  mod_data$resp[,response] <- mod_data$resp[,response] - 0.99999
} 


# Fit models ----------------------------------------------------------------------------------------

# organize data into list to pass to stan
dat_stan <- list(
  n_time  = nrow(mod_data$climate),
  n_lag   = ncol(mod_data$climate),
  y       = mod_data$resp[,response],
  clim    = mod_data$climate,
  clim_means = rowMeans(mod_data$climate),
  clim_yr = list( rowMeans(mod_data$climate[, 1:12]),
                  rowMeans(mod_data$climate[,13:24]),
                  rowMeans(mod_data$climate[,25:36]) ) %>% rbind_l,
  M       = 12,    # number of months in a year
  K       = ncol(mod_data$climate) / 12,
  S       = mod_data$resp$population %>% unique %>% length,
  site_i  = mod_data$resp$population %>% as.factor %>% as.numeric,
  expp_beta = expp_beta
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 3
)

# NULL model (model of the mean)
fit_ctrl1 <- stan(
  file = paste0("stan/",family,"_null.stan"),
  data = dat_stan,
  pars = c('alpha', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# average climate model
fit_ctrl2 <- stan(
  file = paste0("stan/",family,"_ctrl2.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)


# year t
dat_stan$clim_means  <- rowMeans(mod_data$climate[,1:12 ])
fit_yr1 <- stan(
  file = paste0("stan/",family,"_yr.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# year t-1
dat_stan$clim_means  <- rowMeans(mod_data$climate[,13:24])
fit_yr2 <- stan(
  file = paste0("stan/",family,"_yr.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# year t-2
dat_stan$clim_means  <- rowMeans(mod_data$climate[,25:36])
fit_yr3 <- stan(
  file = paste0("stan/",family,"_yr.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)
dat_stan$clim_means  <- rowMeans(mod_data$climate)

# year weights
fit_yr_weight <- stan(
  file = paste0("stan/",family,"_yr_dirichlet.stan"),
  data = dat_stan,
  pars = c('theta', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# year beta (a different beta for each year)
fit_yr_beta <- stan(
  file = paste0("stan/",family,"_yr_beta.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# gaussian moving window
fit_gaus <- stan(
  file = paste0("stan/",family,"_gaus.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# exponential power moving window
fit_expp <- stan(
  file = paste0("stan/",family,"_expp.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# Generalized extreme value 
fit_gev <- stan(
  file = paste0("stan/",family,"_gev.stan"),
  data = dat_stan,
  pars = c('loc', 'scale', "shape", 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# Simplex - 24 months
dat_stan$clim <- t(mod_data$climate)
fit_24 <- stan(
  file = paste0("stan/",family,"_dirichlet.stan"),
  data = dat_stan,
  pars = c('theta', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)
dat_stan$clim <- mod_data$climate
 

# Nested models 
# update data list
dat_stan$clim         <- t(mod_data$climate)
dat_stan$clim1        <- t(mod_data$climate)[1:12 ,]
dat_stan$clim2        <- t(mod_data$climate)[13:24,]
dat_stan$clim3        <- t(mod_data$climate)[25:36,]

# Simplex nested
fit_24_nest <- stan(
  file = paste0("stan/",family,"_dirichlet_nest.stan"),
  data = dat_stan,
  pars = c('theta_y', 'theta_m', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# Generalized Extreme Value nested
fit_gev_nest <- stan(
  file = paste0("stan/",family,"_gev_nest.stan"),
  data = dat_stan,
  pars = c('loc', 'scale', "shape", 'theta_y', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# Power exponential nested 
fit_expp_nest <- stan(
  file = paste0("stan/",family,"_expp_nest.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'theta_y', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)



# parameter values and diagnostics ----------------------------------------------------------------

# list of model fits
mod_fit   <- list(ctrl1   = fit_ctrl1,     ctrl2 = fit_ctrl2,
                  yr1     = fit_yr1,       yr2   = fit_yr2,      yr3   = fit_yr2,
                  yr_wgt  = fit_yr_weight, yr_bet= fit_yr_beta,
                  gaus    = fit_gaus,      expp  = fit_expp,
                  gev     = fit_gev,       simpl = fit_24,
                  simpl_n = fit_24_nest,   gev_n = fit_gev_nest,
                  expp_n  = fit_expp_nest)

# get central tendencies
pars_diag_extract <- function(x){
  
  # central tendencies
  tmp         <- rstan::extract(x) %>% as.data.frame
  tmp         <- setNames(tmp, gsub("\\.","_",names(tmp)))
  par_means   <- sapply(tmp, function(x) mean(x)) %>%
                    setNames( paste0(names(tmp),"_mean") )
  par_medians <- sapply(tmp, function(x) median(x)) %>%
                    setNames( paste0(names(tmp),"_median") )
  central_tend<- c(par_means, par_medians)
  
  # diagnostics
  diverg      <- do.call(rbind, args = get_sampler_params(x, inc_warmup = F))[,5]
  n_diverg    <- length(which(diverg == 1))
  df_summ     <- as.data.frame(summary(x)$summary)
  rhat_high   <- length(which(df_summ$Rhat > 1.1))
  n_eff       <- df_summ$n_eff / length(diverg)
  n_eff_low   <- length(which(n_eff < 0.1))
  mcse_high   <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
  diagnostics <- c(n_diverg = n_diverg, rhat_high = rhat_high,
                   n_eff_low = n_eff_low, mcse_high = mcse_high)
  out         <- c( central_tend, diagnostics ) %>% t %>% as.data.frame
  
  rm(tmp) ; return(out)
  
}

# store posteriors
posterior_extract <- function(model_fit, model_name){
  
  # central tendencies
  tmp     <- rstan::extract(model_fit) %>% as.data.frame
  post_df <- setNames(tmp, gsub("\\.","_",names(tmp)))
  post_df <- tibble::add_column(post_df,
                                model = model_name, .before=1)
  
  rm(tmp) ; return(post_df)
  
}

# calculate central tendencies
pars_diag_l   <- lapply(mod_fit, pars_diag_extract)
mod_pars_diag <- Reduce(function(...) bind_rows(...), pars_diag_l) %>%
                    tibble::add_column(model = names(mod_fit), .before = 1)

# store posteriors
posts_l       <- Map(posterior_extract, mod_fit, names(mod_fit) )
posteriors    <- Reduce(function(...) bind_rows(...), posts_l)


# WAIC model comparison --------------------------------------------------------------------

# wAIC model selection using loo approximation (from library 'loo')
log_liks   <- lapply(mod_fit, extract_log_lik)

# leave-one-out estimates
loo_l      <- lapply(log_liks, loo) %>%
                setNames( c("loo_ctrl1",   "loo_ctrl2", 
                            "loo_yr1",     "loo_yr2",   "loo_yr3", 
                            "loo_yr_wgt",  "loo_yr_bet",
                            "loo_gaus",    "loo_expp",  "loo_gev", "loo_simpl",
                            "loo_simpl_n", "loo_gev_n", "loo_expp_n") )
loo_df     <- loo::compare(loo_l$loo_ctrl1,   loo_l$loo_ctrl2, 
                           loo_l$loo_yr1,     loo_l$loo_yr2,   loo_l$loo_yr3, 
                           loo_l$loo_yr_wgt,  loo_l$loo_yr_bet,   
                           loo_l$loo_gaus,    loo_l$loo_expp,  loo_l$loo_gev,   loo_l$loo_simpl,
                           loo_l$loo_simpl_n, loo_l$loo_gev_n, loo_l$loo_expp_n ) %>%
                as.data.frame %>%
                tibble::add_column(model = gsub("loo_l\\$loo_","",row.names(.) ), .before = 1)

# WAIC estimates
waic_l    <- lapply(log_liks, waic) %>%
                setNames(c("waic_ctrl1",   "waic_ctrl2", 
                           "waic_yr1",     "waic_yr2",   "waic_yr3",  
                           "waic_yr_wgt",  "waic_yr_bet",   
                           "waic_gaus",    "waic_expp",  "waic_gev",   "waic_simpl", 
                           "waic_simpl_n", "waic_gev_n", "waic_expp_n") )
waic_df   <- loo::compare(waic_l$waic_ctrl1,   waic_l$waic_ctrl2, 
                          waic_l$waic_yr1,     waic_l$waic_yr2,     waic_l$waic_yr3,
                          waic_l$waic_yr_wgt,  waic_l$waic_yr_bet,  
                          waic_l$waic_gaus,    waic_l$waic_expp,    waic_l$waic_gev, waic_l$waic_simpl,
                          waic_l$waic_simpl_n, waic_l$waic_gev_n,   waic_l$waic_expp_n) %>%
                as.data.frame %>%
                tibble::add_column(model = gsub("waic_l\\$waic_","",row.names(.) ), .before = 1)


# leave-one-out crossvalidation ------------------------------------------------------------------------

# crossvalidation function
CrossVal <- function(i, mod_data, response) {       # i is index for row to leave out
  
  # identify years
  uniq_yr           <- mod_data$resp$year %>% unique 
  test_i            <- which(mod_data$resp$year == uniq_yr[i])
  
  # put all in matrix form 
  x_clim            <- mod_data$climate
  x_clim_means      <- rowMeans(mod_data$climate)   # climate averages over entire window (for control model #2)
  
  # response variable
  y_train           <- mod_data$resp[-test_i, response]
  y_test            <- mod_data$resp[test_i, response]
  
  # climate variable
  clim_train        <- x_clim[-test_i,]
  clim_test         <- x_clim[test_i,] 
  
  # climate averages over full 24-month window (for control model #2)
  clim_means_train  <- x_clim_means[-test_i]
  clim_means_test   <- x_clim_means[test_i]
  
  # organize data into list to pass to stan
  dat_stan_crossval <- list(
    n_train = length(y_train),  # number of data points in train set (length of response var)
    n_test  = length(y_test),   # number of data points in test set
    n_lag   = ncol(clim_train), # maximum lag
    y_train = array(y_train),
    y_test  = array(y_test),
    clim_train = array(clim_train),
    clim_test  = array(clim_test),
    clim_means_train = array(clim_means_train), # climate averages over full 24-month window (for control model #2)
    clim_means_test  = array(clim_means_test),   # climate averages over full 24-month window (for control model #2)
    clim_yr_train    = list( rowMeans(clim_train[, 1:12]),
                             rowMeans(clim_train[,13:24]),
                             rowMeans(clim_train[,25:36]) ) %>% rbind_l,
    clim_yr_test     = list( rowMeans(clim_test[, 1:12]),
                             rowMeans(clim_test[,13:24]),
                             rowMeans(clim_test[,25:36]) ) %>% rbind_l,
    expp_beta = expp_beta, # beta paramater for exponential power distribution
    M         = 12,  # number of months in a year
    K         = ncol(clim_train) / 12
  )
  
  # fit control 1 (intercept only)
  fit_ctrl1_crossval <- stan(
    file = paste0("stan/",family,"_null_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('alpha', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # fit control 2 (full window climate average)
  fit_ctrl2_crossval <- stan(
    file = paste0("stan/",family,"_ctrl2_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('alpha', 'beta', 'y_sd', 'pred_y', 'log_lik', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # year t
  dat_stan_crossval$clim_means_test   <- rowMeans( mod_data$climate[test_i, 1:12,drop=F] ) %>% array
  dat_stan_crossval$clim_means_train  <- rowMeans( mod_data$climate[-test_i,1:12,drop=F] ) %>% array
  fit_yr1_crossval <- stan(
    file = paste0("stan/",family,"_yr_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('alpha', 'beta', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # year t-1
  dat_stan_crossval$clim_means_test   <- rowMeans( mod_data$climate[test_i, 13:24,drop=F] ) %>% array
  dat_stan_crossval$clim_means_train  <- rowMeans( mod_data$climate[-test_i,13:24,drop=F] ) %>% array
  fit_yr2_crossval <- stan(
    file = paste0("stan/",family,"_yr_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('alpha', 'beta', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
 
  # year t-2
  dat_stan_crossval$clim_means_test   <- rowMeans( mod_data$climate[test_i, 25:36,drop=F] ) %>% array
  dat_stan_crossval$clim_means_train  <- rowMeans( mod_data$climate[-test_i,25:36,drop=F] ) %>% array
  fit_yr3_crossval <- stan(
    file = paste0("stan/",family,"_yr_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('alpha', 'beta', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # year weights
  fit_yr_weight_crossval <- stan(
    file = paste0("stan/",family,"_yr_dirichlet_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('theta', 'alpha', 'beta', 'y_sd', 'log_lik', 'pred_y','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # year beta (a different beta for each year)
  fit_yr_beta_crossval <- stan(
    file = paste0("stan/",family,"_yr_beta_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('alpha', 'beta', 'y_sd', 'log_lik', 'pred_y','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # fit moving window, gaussian
  fit_gaus_crossval <- stan(
    file = paste0("stan/",family,"_gaus_crossval.stan"), 
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains#,
    #control = list(adapt_delta = 0.999, stepsize = 0.005, max_treedepth = 12)
  )
  
  # fit moving window, exponential power
  fit_expp_crossval <- stan(
    file = paste0("stan/",family,"_expp_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains#,
    #control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 12)
  )
  
  # fit GEV (generalized extreme value distribution)
  fit_gev_crossval <- stan(
    file = paste0("stan/",family,"_gev_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('loc','scale','shape','alpha', 'beta', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # update data list
  dat_stan_crossval$clim_train        <- t(dat_stan_crossval$clim_train)
  dat_stan_crossval$clim_test         <- t(dat_stan_crossval$clim_test)
  
  dat_stan_crossval$clim1_train       <- dat_stan_crossval$clim_train[1:12,]
  dat_stan_crossval$clim1_test        <- dat_stan_crossval$clim_test[ 1:12, ,drop=F]
  
  dat_stan_crossval$clim2_train       <- dat_stan_crossval$clim_train[13:24,]
  dat_stan_crossval$clim2_test        <- dat_stan_crossval$clim_test[ 13:24, ,drop=F]

  dat_stan_crossval$clim3_train       <- dat_stan_crossval$clim_train[25:36,]
  dat_stan_crossval$clim3_test        <- dat_stan_crossval$clim_test[ 25:36, ,drop=F]
    
  dat_stan_crossval$clim1_means_train <- colMeans(dat_stan_crossval$clim_train[1:12,])
  dat_stan_crossval$clim1_means_test  <- colMeans(dat_stan_crossval$clim_test[ 1:12, ,drop=F]) %>% array
  
  dat_stan_crossval$clim2_means_train <- colMeans(dat_stan_crossval$clim_train[13:24,])
  dat_stan_crossval$clim2_means_test  <- colMeans(dat_stan_crossval$clim_test[ 13:24, ,drop=F]) %>% array
  
  dat_stan_crossval$clim3_means_train <- colMeans(dat_stan_crossval$clim_train[25:36,])
  dat_stan_crossval$clim3_means_test  <- colMeans(dat_stan_crossval$clim_test[ 25:36, ,drop=F]) %>% array
  
  # fit 24 month simplex
  fit_24_crossval <- stan(
    file = paste0("stan/",family,"_dirichlet_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('theta','alpha', 'beta', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # fit simplex nested within year
  fit_24_nest_crossval <- stan(
    file = paste0("stan/",family,"_dirichlet_nest_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('theta_m', "theta_y",'alpha', 'beta', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # fit gev nested within year
  fit_gev_nest_crossval <- stan(
    file = paste0("stan/",family,"_gev_nest_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('loc','scale','shape', "theta_y",'alpha', 'beta', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # fit exponential power nested within year
  fit_expp_nest_crossval <- stan(
    file = paste0("stan/",family,"_expp_nest_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('sens_mu','sens_sd', "theta_y",'alpha', 'beta', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # posterior mean prediction for the out-of-sample value
  crossval_mods <- list( ctrl1   = fit_ctrl1_crossval, 
                         ctrl2   = fit_ctrl2_crossval, 
                         yr1     = fit_yr1_crossval, 
                         yr2     = fit_yr2_crossval, 
                         yr3     = fit_yr3_crossval, 
                         yr_wgt  = fit_yr_weight_crossval,
                         yr_bet  = fit_yr_beta_crossval,
                         gaus    = fit_gaus_crossval, 
                         expp    = fit_expp_crossval,
                         gev     = fit_gev_crossval,
                         simpl   = fit_24_crossval,
                         simpl_n = fit_24_nest_crossval,
                         gev_n   = fit_gev_nest_crossval,
                         expp_n  = fit_expp_nest_crossval )
                         
  # predictions
  mod_preds <- lapply(crossval_mods, function(x) rstan::extract(x, 'pred_y')$pred_y %>% apply(2,mean) )
  
  # Expected Log Predictive Density
  mod_elpds <- lapply(crossval_mods, function(x){
                                        rstan::extract(x, 'log_lik_test')$log_lik_test %>% 
                                          exp %>%
                                          apply(2,mean) %>%
                                          log 
                                      } )
  
  
  # diagnostics 
  diagnostics <- function(fit_obj, name_mod){
    
    diverg      <- do.call(rbind, args = get_sampler_params(fit_obj, inc_warmup = F))[,5]
    n_diverg    <- length(which(diverg == 1))
    df_summ     <- as.data.frame(summary(fit_obj)$summary)
    rhat_high   <- length(which(df_summ$Rhat > 1.1))
    n_eff       <- df_summ$n_eff / length(diverg)
    n_eff_low   <- length(which(n_eff < 0.1))
    mcse_high   <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
    out         <- data.frame(n_diverg, rhat_high, 
                              n_eff_low, mcse_high) 
    out         <- setNames(out, paste0(names(out),"_",name_mod) )
    return(out)
    
  }
  
  # store diagnostics
  # gaus_expp   <- crossval_mods[c("gaus","expp","gev","simpl")]
  diagnost_l  <- Map(diagnostics, crossval_mods, names(crossval_mods))
  diagnost_df <- do.call(cbind, diagnost_l) %>%
                    bind_cols( unique( dplyr::select(mod_data$resp[test_i,],year) ) )
  
  # function
  pred_elpd_df<- mod_data$resp[test_i,] %>%
                    mutate( # predictions
                            simpl_pred   = mod_preds$simpl,      
                            gev_pred     = mod_preds$gev,
                            gaus_pred    = mod_preds$gaus,
                            expp_pred    = mod_preds$expp,
                            ctrl1_pred   = mod_preds$ctrl1,
                            ctrl2_pred   = mod_preds$ctrl2,
                            yr1_pred     = mod_preds$yr1,
                            yr2_pred     = mod_preds$yr2,
                            yr3_pred     = mod_preds$yr3,
                            yr_wgt_pred  = mod_preds$yr_wgt,
                            yr_bet_pred  = mod_preds$yr_bet,
                            simpl_n_pred = mod_preds$simpl_n,      
                            gev_n_pred   = mod_preds$gev_n,
                            expp_n_pred  = mod_preds$expp_n,
                            
                            # Expected Log Predictive Density
                            simpl_elpd   = mod_elpds$simpl,
                            gev_elpd     = mod_elpds$gev,
                            gaus_elpd    = mod_elpds$gaus,
                            expp_elpd    = mod_elpds$expp,
                            ctrl1_elpd   = mod_elpds$ctrl1,
                            ctrl2_elpd   = mod_elpds$ctrl2,
                            yr1_elpd     = mod_elpds$yr1,
                            yr2_elpd     = mod_elpds$yr2,
                            yr3_elpd     = mod_elpds$yr3,
                            yr_wgt_elpd  = mod_elpds$yr_wgt,
                            yr_bet_elpd  = mod_elpds$yr_bet,
                            simpl_n_elpd = mod_elpds$simpl_n,      
                            gev_n_elpd   = mod_elpds$gev_n,
                            expp_n_elpd  = mod_elpds$expp_n )
                  
  # df to return
  out         <- left_join(pred_elpd_df, diagnost_df)
  
  # remove stanfit objects (garbage collection)
  rm(fit_gev_crossval)        ; rm(fit_24_crossval)
  rm(fit_gaus_crossval)       ; rm(fit_expp_crossval)
  rm(fit_ctrl1_crossval)      ; rm(fit_ctrl2_crossval)
  rm(fit_yr1_crossval)        ; rm(fit_yr2_crossval)   ; rm(fit_yr3_crossval)
  rm(fit_yr_weight_crossval)  ; rm(fit_yr_beta_crossval)
  rm(fit_gev_nest_crossval)
  rm(fit_simpl_nest_crossval)
  rm(fit_expp_nest_crossval)

  return(out)
  
}

# spp-specific cross validation
year_inds   <- seq_along(unique(mod_data$resp$year))
cxval_res   <- lapply( year_inds, CrossVal, mod_data, response)
cxval_pred  <- do.call(rbind, cxval_res) 

# measures of fit -------------------------------------------------------------------------- 

# calculate either mse or deviance
pred_perform <- function(x, mod_data, response, type){
  
  if( type == "mse"){
    res   <- (x - mod_data$resp[,response])^2 %>% mean
  }
  if(type == "deviance"){
    res   <-calc.deviance(x, mod_data$resp[,response], 
                          weights = rep(1, length(x) ),  
                          family="gaussian", calc.mean = TRUE)
  }
  
  return(res)
  
}

# format results into a data frame
perform_format <- function(x, var){
  
  x %>%
    unlist %>%
    t %>% 
    t %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "model") %>%
    mutate( model = gsub("mod_preds.", "", model))  %>%
    setNames( c("model", var) )
  
}

# order of 'mod_data$resp' and 'cxval_pred' need be the same
expect_true( all.equal(dplyr::select(mod_data$resp, year, population),
                       dplyr::select(cxval_pred,    year, population)) )

# mean squared error
mse <- cxval_pred %>% 
          dplyr::select(simpl_pred:expp_n_pred) %>%
          lapply(pred_perform, mod_data, response, "mse") %>%
          perform_format("mse") %>%
          mutate( model = gsub("_pred","",model) )

# # deviance 
# devi <- cxval_pred %>% 
#           dplyr::select(gaus_pred:ctrl2_pred) %>%
#           lapply(pred_perform, mod_data, "deviance") %>%
#           perform_format("deviance") %>%
#           mutate( model = gsub("_pred","",model) )

# Expected Log Predictive Density
elpd <- cxval_pred %>% 
          dplyr::select(simpl_elpd:expp_n_elpd) %>%
          apply(2, sum) %>% 
          as.matrix %>% 
          as.data.frame %>%
          tibble::add_column(model = rownames(.), .before=1) %>%
          mutate( model = gsub("_elpd","",model) ) %>%
          setNames( c("model", "elpd") )

# measures of fit
# mof  <- merge(mse, devi)
mof  <- merge(mse, elpd)

# store results ---------------------------------------------------------------------------
mod_summs <- Reduce(function(...) merge(...), 
                    list(mod_pars_diag, loo_df, waic_df, mof) ) %>%
                    arrange( mse )

write.csv(mod_summs,  paste0(args[3], "_mod_summaries_",spp_name,".csv"), row.names = F)
write.csv(posteriors, paste0(args[3], "_posterior_",spp_name,".csv"), row.names = F)
write.csv(cxval_pred, paste0(args[3], "_crossval_pred_diag_",spp_name,".csv"), row.names = F)
