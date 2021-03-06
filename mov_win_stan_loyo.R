# species-specific lambda and climate data 
setwd("C:/cloud/Dropbox/sAPROPOS project/")
source("~/moving_windows/format_data.R")
library(dplyr)
library(tidyr)
library(loo)
library(rstan)
library(dismo)
library(testthat)

# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# read data -----------------------------------------------------------------------------------------
lam     <- read.csv("DemogData/lambdas_6tr.csv", stringsAsFactors = F) 
m_info  <- read.csv("C:/cloud/MEGA/Projects/sApropos/MatrixEndMonth_information.csv", stringsAsFactors = F)
clim    <- read.csv("DemogData/precip_fc_demo.csv",  stringsAsFactors = F) #%>%
              #mutate( ppt = ppt / 30)
spp     <- clim$species %>% unique

# add monthg info to lambda information
month_add <- m_info %>%
              mutate(SpeciesAuthor = trimws(SpeciesAuthor) ) %>%
              dplyr::select(SpeciesAuthor, MatrixEndMonth)
lam_add   <- subset(lam, SpeciesAuthor %in% m_info$SpeciesAuthor) %>%  
              dplyr::select(-MatrixEndMonth) %>%
              inner_join(month_add)
lam_min   <- subset(lam, SpeciesAuthor %in% setdiff(lam$SpeciesAuthor, m_info$SpeciesAuthor) )
lambdas   <- bind_rows(lam_min, lam_add) %>%
              subset( !is.na(MatrixEndMonth) )


# format data ---------------------------------------------------------------------------------------
spp       <- lambdas$SpeciesAuthor %>% unique %>% .[-3]   
for(ii in 1:35){
spp_name      <- spp[ii] # test run w/ spp number 1
m_back        <- 24     # months back
expp_beta     <- 20

# lambda data
spp_lambdas   <- format_species(spp_name, lambdas)

# climate data
clim_separate <- clim_list(spp_name, clim, spp_lambdas)
clim_detrnded <- lapply(clim_separate, clim_detrend)
clim_mats     <- Map(clim_long, clim_detrnded, spp_lambdas, m_back)

# model data
mod_data          <- lambda_plus_clim(spp_lambdas, clim_mats)
mod_data$climate  <- mod_data$climate #/ diff(range(mod_data$climate))


# model fits ---------------------------------------------------------------------------------

# organize data into list to pass to stan
dat_stan <- list(
  n_time = nrow(mod_data$climate),
  n_lag  = ncol(mod_data$climate),
  y      = mod_data$lambdas$log_lambda,
  clim   = mod_data$climate ,
  clim_means = rowMeans(mod_data$climate),
  expp_beta = expp_beta
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 4
)

# moving window, gaussian
fit_gaus <- stan( 
  file = 'Code/climwin/stan_movwin/stan/movwin_gaus.stan',
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# moving window, exponential power
fit_expp <- stan( 
  file = 'Code/climwin/stan_movwin/stan/movwin_expp.stan',
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# control 1 (intercept only)
fit_ctrl1 <- stan(
  file = 'Code/climwin/stan_movwin/stan/movwin_ctrl1.stan',
  data = dat_stan,
  pars = c('alpha', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)

# control 2 (full window average)
fit_ctrl2 <- stan(
  file = 'Code/climwin/stan_movwin/stan/movwin_ctrl2.stan',
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)


# parameter values and diagnostics ----------------------------------------------------------------

# list of model fits
mod_fit   <- list(gaus = fit_gaus, expp = fit_expp, 
                  ctrl1 = fit_ctrl1, ctrl2 = fit_ctrl2)

# parameter values
pars      <- c('sens_mu', 'sens_sd', 'alpha', 'beta', "y_sd")

# get central tendencies
pars_diag_extract <- function(x){
  
  # central tendencies
  tmp         <- rstan::extract(x)
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
  tmp         <- rstan::extract(model_fit)
  post_df     <- do.call(cbind, tmp) %>% as.data.frame
  ll_id       <- grep("V", colnames(post_df) )
  new_names   <- paste0("log_lik_", 1:length(ll_id) )
  names(post_df)[ll_id] <- new_names # no way to do this in dplyr
  post_df     <- tibble::add_column(post_df, 
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
log_liks  <- lapply(mod_fit, extract_log_lik)
# shiny_eval<- lapply(mod_fit[[1]], launch_shinystan) # evaluate model convergence and fit using library 'shinystan'

# leave-one-out estimates
loo_l      <- lapply(log_liks, loo) %>%
                setNames(c("loo_gaus", "loo_expp", "loo_ctrl1",  "loo_ctrl2"))
loo_df     <- loo::compare(loo_l$loo_gaus, loo_l$loo_expp, loo_l$loo_ctrl1, loo_l$loo_ctrl2) %>%
                as.data.frame %>%
                tibble::add_column(model = gsub("loo_","",names(loo_l) ), .before = 1)

# WAIC estimates
waic_l    <- lapply(log_liks, waic) %>%
                setNames(c("waic_gaus", "waic_expp", "waic_ctrl1",  "waic_ctrl2"))
waic_df   <- loo::compare(waic_l$waic_gaus, waic_l$waic_expp, waic_l$waic_ctrl1, waic_l$waic_ctrl2) %>%
                as.data.frame %>%
                tibble::add_column(model = gsub("waic_","",names(waic_l) ), .before = 1)

# leave-one-YEAR-out crossvalidation ---------------------------------------------------

# crossvalidation function
CrossVal <- function(i, mod_data) {       # i is index for row to leave out
  
  # identify years
  uniq_yr           <- mod_data$lambdas$year %>% unique 
  test_i            <- which(mod_data$lambdas$year == uniq_yr[i])
  
  # put all in matrix form 
  x_clim            <- mod_data$climate
  x_clim_means      <- rowMeans(mod_data$climate)   # climate averages over entire window (for control model #2)
  
  # response variable
  y_train           <- mod_data$lambdas$log_lambda[-test_i]
  y_test            <- mod_data$lambdas$log_lambda[test_i]
  
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
    clim_test = array(clim_test),
    clim_means_train = array(clim_means_train), # climate averages over full 24-month window (for control model #2)
    clim_means_test = array(clim_means_test),   # climate averages over full 24-month window (for control model #2)
    expp_beta = expp_beta       # beta paramater for exponential power distribution
  )
  
  # fit moving window, gaussian
  fit_gaus_crossval <- stan( 
    file = 'Code/climwin/stan_movwin/stan/movwin_gaus_loyo.stan',
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik', 'pred_y'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.999, stepsize = 0.005, max_treedepth = 12)
  )
  
  # fit moving window, exponential power
  fit_expp_crossval <- stan( 
    file = 'Code/climwin/stan_movwin/stan/movwin_expp_loyo.stan',
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik', 'pred_y'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 12)
  )
  
  # fit control 1 (intercept only)
  fit_ctrl1_crossval <- stan(
    file = 'Code/climwin/stan_movwin/stan/movwin_ctrl1_loyo.stan',
    data = dat_stan_crossval,
    pars = c('alpha', 'y_sd', 'log_lik', 'pred_y'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # fit control 2 (full window climate average)
  fit_ctrl2_crossval <- stan(
    file = 'Code/climwin/stan_movwin/stan/movwin_ctrl2_loyo.stan',
    data = dat_stan_crossval,
    pars = c('alpha', 'beta', 'y_sd', 'log_lik', 'pred_y'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # posterior mean prediction for the out-of-sample value
  crossval_mods <- list(gaus = fit_gaus_crossval, 
                        expp = fit_expp_crossval, 
                        ctrl1 = fit_ctrl1_crossval, 
                        ctrl2 = fit_ctrl2_crossval
                        )
  
  # predictions
  mod_preds <- lapply(crossval_mods, function(x) rstan::extract(x, 'pred_y')$pred_y %>% apply(2,mean) )

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
  gaus_expp   <- crossval_mods[1:2]
  diagnost_l  <- Map(diagnostics, gaus_expp, names(gaus_expp))
  diagnost_df <- do.call(cbind, diagnost_l) %>%
                    bind_cols( unique( dplyr::select(mod_data$lambdas[test_i,],year) ) )
  
  # function
  pred_df     <- mod_data$lambdas[test_i,] %>%
                    mutate( gaus_pred  = mod_preds$gaus,
                            expp_pred  = mod_preds$expp,
                            cltr1_pred = mod_preds$ctrl1,
                            cltr2_pred = mod_preds$ctrl2 )
  
  # df to return
  out         <- left_join(pred_df, diagnost_df)
         
  # remove stanfit objects (garbage collection)
  rm(fit_gaus_crossval)
  rm(fit_expp_crossval)
  rm(fit_ctrl1_crossval)
  rm(fit_ctrl2_crossval)
  
  return(out)
  
}

# spp-specific cross validation
year_inds   <- seq_along(unique(mod_data$lambdas$year))
cxval_res   <- lapply( year_inds, CrossVal, mod_data)
cxval_pred  <- do.call(rbind, cxval_res) 


# measures of fit -------------------------------------------------------------------------- 

# calculate either mse or deviance
pred_perform <- function(x, mod_data, type){
  
  if( type == "mse"){
    res   <- (x - mod_data$lambdas$log_lambda)^2 %>% mean
  }
  if(type == "deviance"){
    res   <-calc.deviance(mod_data$lambdas$log_lambda, x, 
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

# mean squared error
mse <- cxval_pred %>% 
          dplyr::select( gaus_pred:cltr2_pred ) %>%
          lapply(pred_perform, mod_data, "mse") %>%
          perform_format("mse")

# deviance 
devi <- cxval_pred %>% 
          dplyr::select( gaus_pred:cltr2_pred ) %>%
          lapply(pred_perform, mod_data, "deviance") %>%
          perform_format("deviance")

# measures of fit
mof  <- merge(mse, devi)


# store results ---------------------------------------------------------------------------
mod_summs <- Reduce(function(...) merge(...), 
                    list(mod_pars_diag, loo_df, waic_df, mof) ) %>%
                    arrange( mse )
write.csv(mod_summs, paste0("C:/cloud/MEGA/Projects/sApropos/results/mod_summaries_",spp_name,".csv"), row.names = F)
write.csv(posteriors, paste0("C:/cloud/MEGA/Projects/sApropos/results/posterior_",spp_name,".csv"), row.names = F)
write.csv(cxval_pred, paste0("C:/cloud/MEGA/Projects/sApropos/results/crossval_pred_diag_",spp_name,".csv"), row.names = F)
}

