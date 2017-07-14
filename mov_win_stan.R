# species-specific lambda and climate data 
setwd("C:/cloud/Dropbox/sAPROPOS project/")
source("~/moving_windows/format_data.R")
library(dplyr)
library(tidyr)
library(loo)
library(rstan)
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
              select(SpeciesAuthor, MatrixEndMonth)
lam_add   <- subset(lam, SpeciesAuthor %in% m_info$SpeciesAuthor) %>%  
              select(-MatrixEndMonth) %>%
              inner_join(month_add)
lam_min   <- subset(lam, SpeciesAuthor %in% setdiff(lam$SpeciesAuthor, m_info$SpeciesAuthor) )
lambdas   <- bind_rows(lam_min, lam_add) %>%
              subset( !is.na(MatrixEndMonth) )


# format data ---------------------------------------------------------------------------------------
spp_name      <- spp[4] # test run w/ spp number 1
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
  control = list(adapt_delta = 0.999, stepsize = 0.005, max_treedepth = 12)
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
  control = list(adapt_delta = 0.999, stepsize = 0.005, max_treedepth = 12)
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


# store parameter values ----------------------------------------------------------------

# list of model fits
mod_fit   <- list(gaus = fit_gaus, expp = fit_expp, 
                  ctrl1 = fit_ctrl1, ctrl2 = fit_ctrl2)

# parameter values
pars      <- c('sens_mu', 'sens_sd', 'alpha', 'beta', "y_sd")

# get central tendencies
central_tend_get <- function(x){
  
  tmp         <- rstan::extract(x)
  par_means   <- sapply(tmp, function(x) mean(x)) %>% 
                    setNames( paste0(names(tmp),"_mean") ) 
  par_medians <- sapply(tmp, function(x) median(x)) %>%
                    setNames( paste0(names(tmp),"_median") )
  out         <- c(par_means, par_medians) %>% t %>% as.data.frame
  
  rm(tmp) ; return(out)
}

# calculate and store central tendencies
centr_tend_list     <- lapply(mod_fit, central_tend_get)
central_tendencies  <- Reduce(function(...) bind_rows(...), centr_tend_list) %>%
                          tibble::add_column(model = names(mod_fit), .before = 1)
#write.csv(central_tendencies, paste0(spp_name,"_central_tendencies.csv"), row.names=F)

  
# WAIC model comparison --------------------------------------------------------------------

# # Rhat convergence check  
# rhat_check <- function(input_mod){
#   
#   rhats <- summary(input_mod)$summary[,"Rhat"]
#   if( any(rhats > 1.1) ){
#     TRUE
#   } else{ FALSE }
#   
# }
# expect_false( any(sapply(mod_fit, rhat_check)) )


# wAIC model selection using loo approximation (from library 'loo')
log_liks  <- lapply(mod_fit, extract_log_lik)
loos      <- lapply(log_liks, loo) %>%
              setNames(c("loo_gaus", "loo_expp", "loo_ctrl1",  "loo_ctrl2"))
# shiny_eval <- lapply(fits, launch_shinystan) # evaluate model convergence and fit using library 'shinystan'
waics     <- compare(loos$loo_gaus, loos$loo_expp, loos$loo_ctrl1, loos$loo_ctrl2)
#write.csv(waics, paste0("waic_",spp_name,".csv"), row.names=F)




# leave-one-out crossvalidation ---------------------------------------------------

# crossvalidation function
CrossVal <- function(i, mod_data) {       # i is index for row to leave out
  
  # put all in matrix form 
  x_clim            <- mod_data$climate
  x_clim_means      <- rowMeans(mod_data$climate)   # climate averages over entire window (for control model #2)
  
  # response variable
  y_train           <- mod_data$lambdas$log_lambda[-i]
  y_test            <- mod_data$lambdas$log_lambda[i]
  
  # climate variable
  clim_train        <- x_clim[-i,]
  clim_test         <- x_clim[i,] %>% as.numeric
  
  # climate averages over full 24-month window (for control model #2)
  clim_means_train  <- x_clim_means[-i]
  clim_means_test   <- x_clim_means[i] 
  
  # organize data into list to pass to stan
  dat_stan_crossval <- list(
    n_time = nrow(clim_train),  # number of years (length of response var)
    n_lag = ncol(clim_train),   # maximum lag
    y_train = y_train,
    y_test = y_test,
    clim_train = clim_train,
    clim_test = clim_test,
    clim_means_train = clim_means_train,           # climate averages over full 24-month window (for control model #2)
    clim_means_test = clim_means_test,            # climate averages over full 24-month window (for control model #2)
    expp_beta = expp_beta       # beta paramater for exponential power distribution
  )
  
  # fit moving window, gaussian
  fit_gaus_crossval <- stan( 
    file = 'Code/climwin/stan_movwin/stan/movwin_gaus_crossval.stan',
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik', 'pred_y'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.995, stepsize = 0.005, max_treedepth = 12)
  )
  
  # fit moving window, exponential power
  fit_expp_crossval <- stan( 
    file = 'Code/climwin/stan_movwin/stan/movwin_expp_crossval.stan',
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik', 'pred_y'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains,
    control = list(adapt_delta = 0.995, stepsize = 0.1, max_treedepth = 12)
  )
  
  # fit control 1 (intercept only)
  fit_ctrl1_crossval <- stan(
    file = 'Code/climwin/stan_movwin/stan/movwin_ctrl1_crossval.stan',
    data = dat_stan_crossval,
    pars = c('alpha', 'y_sd', 'log_lik', 'pred_y'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # fit control 2 (full window climate average)
  fit_ctrl2_crossval <- stan(
    file = 'Code/climwin/stan_movwin/stan/movwin_ctrl2_crossval.stan',
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
  mod_preds <- lapply(crossval_mods, function(x) rstan::extract(x, 'pred_y')$pred_y %>% mean )

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
    
  }
  
  # store diagnostics
  gaus_expp   <- crossval_mods[1:2]
  diagnost_l  <- Map(diagnostics, gaus_expp, names(gaus_expp))
  diagnost_df <- do.call(cbind, diagnost_l)        
  
  # df to return
  out <- data.frame(mod_preds$gaus, mod_preds$expp, mod_preds$ctrl1, mod_preds$ctrl2) %>%
            bind_cols(diagnost_df)
                    
  # remove stanfit objects (garbage collection)
  rm(fit_gaus_crossval)
  rm(fit_expp_crossval)
  rm(fit_ctrl1_crossval)
  rm(fit_ctrl2_crossval)
  
  return(out)
  
}

# spp-specific cross validation
cxval_res   <- lapply(1:nrow(mod_data$lambdas), CrossVal, mod_data) %>%
                rbindlist(idcol=T)
cxval_pred  <- do.call(rbind, cxval_res) 

# mean squared errors
mses         <- cxval_pred %>% 
                  select(mod_preds.gaus:mod_preds.ctrl2) %>%
                  sweep(1, mod_data$lambdas$log_lambda, "-") %>%
                  .^2 %>%
                  colMeans
