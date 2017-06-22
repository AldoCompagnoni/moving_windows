
##### preliminaries ---------------------------------------------------

# load libraries
library(dplyr)
library(tibble)
library(tidyr)
library(data.table)
library(rstan)
library(shinystan)
library(loo)
library(ggplot2)
library(gridExtra)
library(caret)

# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# setwd
setwd('~/desktop/stan_movwin/') # ***change to relevant path***

# read Dalgleish et al. data
d <- read.csv("data/lambdas.csv")
precip <- read.csv("data/monthly_ppt_Dalgleish.csv")

# analysis parameters  
spp <- 1
m_back <- 24
expp_beta <- 20  # beta param for exponential power function (higher means more square)

# exponential power distribution (aka generalized normal)
dexppow <- function(x, mu, sigma, beta) {
  return((beta / (2 * sigma * gamma(1/beta)) ) * exp(-(abs(x - mu)/sigma)^beta))
}




##### format species ---------------------------------------------------

# select species
spp_dur <- d %>% 
  group_by(SpeciesAccepted) %>% 
  summarise(duration = length(unique(MatrixStartYear)))

# species list
spp_list <- d %>%
  subset(Lat == 38.8 & Lon == -99.2) %>%
  .[,"SpeciesAccepted"] %>%
  as.character %>%
  unique

# lambdas, all species
xx_full <- d %>%
  subset(SpeciesAccepted %in% spp_list) %>%
  select(SpeciesAccepted, MatrixEndYear, lambda) %>%
  setNames(c("species", "year","lambda")) %>%
  mutate(log_lambda = log(lambda), species = as.character(species))

# lambdas, focal species
xx <- xx_full %>% subset(species == spp_list[spp])




##### format climate ---------------------------------------------------------

# detrend climate
m_means <- colMeans(precip, na.rm=T)[-1]
d_precip <- apply(precip[,-1], 2, FUN = scale, center = T, scale = T) 
det_precip <- cbind(precip[,1], d_precip) %>% as.data.frame
names(det_precip)[1] <- "YEAR"           

# select precipitation range
precip_long <- det_precip %>%
  subset(YEAR < 1974 & YEAR > 1934) %>%
  gather(month, precip, JAN:DEC) %>%
  setNames(c("year", "month", "precip")) %>% 
  mutate(month_num = factor(month, levels = toupper(month.abb))) %>% 
  mutate(month_num = as.numeric(month_num)) %>% 
  arrange(year, month_num)

# array for number of months before each sampling date
precip_form <- function(x, dat, var) {
  id <- which(dat$year == x & dat$month == "JUL")
  r  <- c(id:(id - (m_back-1)))
  return(dat[r,var])
}

# calculate monthly precipitation values
years   <- xx$year %>% unique() %>% sort()
precip_l  <- lapply(years, precip_form, precip_long, "precip")

# put all in matrix form 
mat_form<- function(x, years) { 
  do.call(cbind, x) %>% 
    as.data.frame %>%
    setNames(years) 
}

x_precip  <- t(mat_form(precip_l, years))
x_precip_means <- rowMeans(x_precip)      # climate averages over entire window (for control model #2)



##### model ---------------------------------------------------------

# organize data into list to pass to stan
dat_stan <- list(
  n_time = nrow(x_precip),
  n_lag = ncol(x_precip),
  y = xx$log_lambda,
  clim = x_precip,
  clim_means = x_precip_means,
  expp_beta = expp_beta
)

# moving window, gaussian
fit_gaus <- stan( 
  file = 'stan/movwin_gaus.stan',
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = 1000,
  iter = 3000,
  thin = 2,
  chains = 2,
  control = list(adapt_delta = 0.995, stepsize = 0.005, max_treedepth = 12)
)

# moving window, exponential power
fit_expp <- stan( 
  file = 'stan/movwin_expp.stan',
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = 1000,
  iter = 3000,
  thin = 2,
  chains = 2,
  control = list(adapt_delta = 0.9, stepsize = 0.1, max_treedepth = 12)
)

# control 1 (intercept only)
fit_ctrl1 <- stan(
  file = 'stan/movwin_ctrl1.stan',
  data = dat_stan,
  pars = c('alpha', 'y_sd', 'log_lik'),
  warmup = 1000,
  iter = 3000,
  thin = 2,
  chains = 2
)

# control 2 (full window average)
fit_ctrl2 <- stan(
  file = 'stan/movwin_ctrl2.stan',
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = 1000,
  iter = 3000,
  thin = 2,
  chains = 2
)

# evaluate model convergence and fit using library 'shinystan'
# launch_shinystan(fit_guas)
# launch_shinystan(fit_expp)
# launch_shinystan(fit_ctrl1)
# launch_shinystan(fit_ctrl2)

# evaluate models with loo approximation (from library 'loo')
ll_gaus <- extract_log_lik(fit_gaus)
ll_expp <- extract_log_lik(fit_expp)
ll_ctrl1 <- extract_log_lik(fit_ctrl1)
ll_ctrl2 <- extract_log_lik(fit_ctrl2)

loo_gaus <- loo(ll_gaus)
loo_expp <- loo(ll_expp)
loo_ctrl1 <- loo(ll_ctrl1)
loo_ctrl2 <- loo(ll_ctrl2)

compare(loo_gaus, loo_expp, loo_ctrl1, loo_ctrl2)

# extract posterior parameter estimates from fitted stan models
sens_mu_gaus <- rstan::extract(fit_gaus, 'sens_mu')$sens_mu
sens_sd_gaus <- rstan::extract(fit_gaus, 'sens_sd')$sens_sd
alpha_gaus <- rstan::extract(fit_gaus, 'alpha')$alpha
beta_gaus <- rstan::extract(fit_gaus, 'beta')$beta

sens_mu_expp <- rstan::extract(fit_expp, 'sens_mu')$sens_mu
sens_sd_expp <- rstan::extract(fit_expp, 'sens_sd')$sens_sd
alpha_expp <- rstan::extract(fit_expp, 'alpha')$alpha
beta_expp <- rstan::extract(fit_expp, 'beta')$beta

# posterior medians
med_sens_mu_gaus <- median(sens_mu_gaus)
med_sens_sd_gaus <- median(sens_sd_gaus)
med_alpha_gaus <- median(alpha_gaus)
med_beta_gaus <- median(beta_gaus)

med_sens_mu_expp <- median(sens_mu_expp)
med_sens_sd_expp <- median(sens_sd_expp)
med_alpha_expp <- median(alpha_expp)
med_beta_expp <- median(beta_expp)

# predicted sensitivity functions
t_sens <- seq(1, m_back, 1)
sens_gaus <- dnorm(t_sens, med_sens_mu_gaus, med_sens_sd_gaus)
sens_gaus <- sens_gaus / sum(sens_gaus)
sens_expp <- sapply(t_sens, FUN = dexppow, mu = med_sens_mu_expp, sigma = med_sens_sd_gaus, beta = expp_beta)
sens_expp <- sens_expp / sum(sens_expp)
sens_df <- data.frame(time_rel = 1:m_back, sens_gaus, sens_expp)

# predicted sensitivity-weighted precipitation values
precip_sens_df <- x_precip %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'year') %>%
  gather_('time_rel', 'precip', paste0('V', 1:m_back)) %>% 
  arrange(year) %>% 
  mutate(time_rel = as.numeric(gsub('V', '', time_rel))) %>% 
  left_join(sens_df, by = 'time_rel') %>% 
  group_by(year) %>% 
  summarize(precip_gaus = sum(sens_gaus * precip), precip_expp = sum(sens_expp * precip)) %>% 
  mutate(year = as.integer(year)) %>% 
  left_join(xx, by = 'year')


# generate prediction lines for each posterior sample (for sensitivity functions)
t_sens_plot <- seq(1, m_back, length.out = 250)  # generate evenly-spaced time values over relevant range
sens_full_gaus <- mapply(function(x, y) dnorm(t_sens_plot, x, y), x = sens_mu_expp, y = sens_sd_expp)
sens_full_expp <- mapply(function(x, y) sapply(t_sens_plot, FUN = dexppow, x, y, expp_beta), x = sens_mu_expp, y = sens_sd_expp)


# generate prediction lines for each posterior sample (for relationship between response and climate)
precip_pred <- pretty(c(precip_sens_df$precip_gaus, precip_sens_df$precip_expp), 40)  # generate evenly-spaced precipitation values over relevant range
beta_full_gaus <- mapply(function(a, b) a + b * precip_pred, a = alpha_gaus, b = beta_gaus)
beta_full_expp <- mapply(function(a, b) a + b * precip_pred, a = alpha_expp, b = beta_expp)

# arrange prediction lines for each posterior sample in tidy data frame
sens_df <- sens_full_gaus %>% 
  as.data.frame() %>% 
  mutate(t_sens_plot = t_sens_plot) %>% 
  melt(id.vars = 't_sens_plot', value.name = 'val_gaus') %>% 
  mutate(val_expp = c(sens_full_expp)) %>% 
  as_tibble()

# arrange prediction lines for each posterior sample in tidy data frame
sens_df <- sens_full_gaus %>% 
  as.data.frame() %>% 
  mutate(t_sens_plot = t_sens_plot) %>% 
  melt(id.vars = 't_sens_plot', value.name = 'val_gaus') %>% 
  mutate(val_expp = c(sens_full_expp)) %>% 
  as_tibble()

pred_df <- beta_full_gaus %>% 
  as.data.frame() %>% 
  mutate(precip_pred = precip_pred) %>% 
  melt(id.vars = 'precip_pred', value.name = 'val_gaus') %>% 
  mutate(val_expp = c(beta_full_expp)) %>% 
  as_tibble()

# generate best-fit prediction lines
sens_best_df <- data.frame(
  t_sens_plot,
  sens_best_gaus = dnorm(t_sens_plot, med_sens_mu_gaus, med_sens_sd_gaus),
  sens_best_expp = sapply(t_sens_plot, FUN = dexppow, mu = med_sens_mu_expp, sigma = med_sens_sd_gaus, beta = expp_beta)
) # note that we don't need to scale sensitivity values by their sum for graphing purposes

pred_best_df <- data.frame(
  precip_pred,
  y_pred_gaus = med_alpha_gaus + med_beta_gaus * precip_pred,
  y_pred_expp = med_alpha_expp + med_beta_expp * precip_pred
)

# create plots
tt <- theme_bw() +    # custom ggplot theme
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 13),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(.3, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .3, 0, 0, unit = 'cm')))

p1 <- ggplot(sens_best_df) +
  geom_line(aes(t_sens_plot-1, sens_best_gaus), col = 'red', size = 1.5) +
  ggtitle('Gaussian') +
  xlab(NULL) + ylab(NULL) +
  tt + theme(axis.text = element_blank(), axis.ticks.y = element_blank())

p2 <- ggplot(sens_best_df) +
  geom_line(aes(t_sens_plot-1, sens_best_expp), col = 'red', size = 1.5) +
  ggtitle('Exponential Power') +
  xlab(NULL) + ylab(NULL) +
  tt + theme(axis.text = element_blank(), axis.ticks.y = element_blank())

p3 <- ggplot(sens_df, aes(t_sens_plot-1, val_gaus, group = variable)) +
  geom_line(alpha = 0.06) +
  xlab('Months prior to demographic sampling') + ylab('Climate sensitivity') +
  tt

p4 <- ggplot(sens_df, aes(t_sens_plot-1, val_expp, group = variable)) +
  geom_line(alpha = 0.06) +
  xlab('Months prior to demographic sampling') + ylab('Climate sensitivity') +
  tt

p5 <- ggplot(pred_df, aes(precip_pred, val_gaus, group = variable)) +
  geom_line(alpha = 0.04) +
  geom_line(data = pred_best_df, inherit.aes = F, aes(precip_pred, y_pred_gaus), col = 'red', size = 2) +
  geom_point(data = precip_sens_df, inherit.aes = F, aes(precip_gaus, log_lambda), col = 'blue', size = 4, alpha = 0.5) +
  xlab('Weighted precipitation') + ylab('ln Lambda') +
  tt

p6 <- ggplot(pred_df, aes(precip_pred, val_expp, group = variable)) +
  geom_line(alpha = 0.04) +
  geom_line(data = pred_best_df, inherit.aes = F, aes(precip_pred, y_pred_expp), col = 'red', size = 2) +
  geom_point(data = precip_sens_df, inherit.aes = F, aes(precip_expp, log_lambda), col = 'blue', size = 4, alpha = 0.5) +
  xlab('Weighted precipitation') + ylab('ln Lambda') +
  tt

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g4 <- ggplotGrob(p4)
g5 <- ggplotGrob(p5)
g6 <- ggplotGrob(p6)

g1$widths <- g3$widths
g5$widths <- g3$widths
g2$widths <- g4$widths
g6$widths <- g4$widths

mw_plot <- arrangeGrob(g1, g2, g3, g4, g5, g6, nrow = 3, heights = c(0.2, 0.4, 0.4))

# plot to device
dev.off()
quartz(height = 10, width = 15)
grid.arrange(mw_plot)

# save image to file
file_out <- paste0('analysis/stan_movwin_spp_', spp, '.tiff')
ggsave(file_out, mw_plot, height = 10, width = 15, units = 'in', dpi = 300)




##### cross validation ---------------------------------------------------------

CrossVal <- function(species_focal, i) {       # i is index for row to leave out
  
  # lambdas, focal species
  xx <- xx_full %>% subset(species == species_focal)
  
  # calculate monthly precipitation values
  years   <- xx$year %>% unique() %>% sort()
  precip_l  <- lapply(years, precip_form, precip_long, "precip")
  
  # put all in matrix form 
  x_precip  <- t(mat_form(precip_l, years))
  x_precip_means <- rowMeans(x_precip)   # climate averages over entire window (for control model #2)
  
  # response variable
  y_train <- xx$log_lambda[-i]
  y_test <- xx$log_lambda[i]
  
  # climate variable
  clim_train  <- x_precip[-i,]
  clim_test  <- x_precip[i,]
  
  # climate averages over full 24-month window (for control model #2)
  clim_means_train <- x_precip_means[-i]
  clim_means_test <- x_precip_means[i]
  
  # organize data into list to pass to stan
  dat_stan_crossval <- list(
    n_time = nrow(clim_train),  # number of years (length of response var)
    n_lag = ncol(clim_train),   # maximum lag
    y_train,
    y_test,
    clim_train,
    clim_test,
    clim_means_train,           # climate averages over full 24-month window (for control model #2)
    clim_means_test,            # climate averages over full 24-month window (for control model #2)
    expp_beta = expp_beta       # beta paramater for exponential power distribution
  )
  
  # fit moving window, gaussian
  fit_gaus_crossval <- stan( 
    file = 'stan/movwin_gaus_crossval.stan',
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik', 'pred_y'),
    warmup = 1000,
    iter = 3000,
    thin = 2,
    chains = 2,
    control = list(adapt_delta = 0.995, stepsize = 0.005, max_treedepth = 12)
  )
  
  # fit moving window, exponential power
  fit_expp_crossval <- stan( 
    file = 'stan/movwin_expp_crossval.stan',
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik', 'pred_y'),
    warmup = 1000,
    iter = 3000,
    thin = 2,
    chains = 2,
    control = list(adapt_delta = 0.9, stepsize = 0.1, max_treedepth = 12)
  )
  
  # fit control 1 (intercept only)
  fit_ctrl1_crossval <- stan(
    file = 'stan/movwin_ctrl1_crossval.stan',
    data = dat_stan_crossval,
    pars = c('alpha', 'y_sd', 'log_lik', 'pred_y'),
    warmup = 1000,
    iter = 3000,
    thin = 2,
    chains = 2
  )
  
  # fit control 2 (full window climate average)
  fit_ctrl2_crossval <- stan(
    file = 'stan/movwin_ctrl2_crossval.stan',
    data = dat_stan_crossval,
    pars = c('alpha', 'beta', 'y_sd', 'log_lik', 'pred_y'),
    warmup = 1000,
    iter = 3000,
    thin = 2,
    chains = 2
  )
  
  # posterior mean prediction for the out-of-sample value
  pred_gaus <- mean(rstan::extract(fit_gaus_crossval, 'pred_y')$pred_y)
  pred_expp <- mean(rstan::extract(fit_expp_crossval, 'pred_y')$pred_y)
  pred_c1 <- mean(rstan::extract(fit_ctrl1_crossval, 'pred_y')$pred_y)
  pred_c2 <- mean(rstan::extract(fit_ctrl2_crossval, 'pred_y')$pred_y)
  
  # diagnostics for gaussian moving window
  diverg_gaus <- do.call(rbind, args = get_sampler_params(fit_gaus_crossval, inc_warmup = F))[,5]
  n_diverg_gaus <- length(which(diverg_gaus == 1))
  df_summ_gaus <- as.data.frame(summary(fit_gaus_crossval)$summary)
  rhat_high_gaus <- length(which(df_summ_gaus$Rhat > 1.1))
  n_eff_gaus <- df_summ_gaus$n_eff / length(diverg_gaus)
  n_eff_low_gaus <- length(which(n_eff_gaus < 0.1))
  mcse_high_gaus <- length(which(df_summ_gaus$se_mean / df_summ_gaus$sd > 0.1))
  
  # diagnostics for exponential-power moving window
  diverg_expp <- do.call(rbind, args = get_sampler_params(fit_expp_crossval, inc_warmup = F))[,5]
  n_diverg_expp <- length(which(diverg_expp == 1))
  df_summ_expp <- as.data.frame(summary(fit_expp_crossval)$summary)
  rhat_high_expp <- length(which(df_summ_expp$Rhat > 1.1))
  n_eff_expp <- df_summ_expp$n_eff / length(diverg_expp)
  n_eff_low_expp <- length(which(n_eff_expp < 0.1))
  mcse_high_expp <- length(which(df_summ_expp$se_mean / df_summ_expp$sd > 0.1))
  
  # df to return
  out <- data.frame(pred_gaus, pred_expp, pred_c1, pred_c2,
                    n_diverg_gaus, rhat_high_gaus, n_eff_low_gaus, mcse_high_gaus,
                    n_diverg_expp, rhat_high_expp, n_eff_low_expp, mcse_high_expp)
  
  # remove stanfit objects (garbage collection)
  rm(fit_gaus_crossval)
  rm(fit_expp_crossval)
  rm(fit_ctrl1_crossval)
  rm(fit_ctrl2_crossval)
  
  return(out)
}


# # perform leave-one-out cross validation for all 5 species and all 4 models (takes about 1hr on my Macbook)
# crossval_results <- xx_full %>%
#   group_by(species) %>%
#   mutate(index = 1:n()) %>%
#   ungroup() %>%
#   group_by(species, index, year) %>%
#   do(CrossVal(species_focal = .$species, i = .$index)) %>% ungroup()
# 
# 
# # write results to file
# write.csv(crossval_results, 'analysis/stan_movwin_crossval.csv', row.names = F)


# read crossval results for each species
crossval_results <- read.csv('analysis/stan_movwin_crossval.csv', stringsAsFactors = F)

# summarize crossval results
crossval_summary <- crossval_results %>% 
  left_join(dplyr::select(xx_full, species, year, log_lambda), by = c('species', 'year')) %>% 
  group_by(species) %>% 
  summarize(rmse_guas = sqrt(mean((pred_gaus - log_lambda)^2)),
            rmse_expp = sqrt(mean((pred_expp - log_lambda)^2)),
            rmse_ctrl1 = sqrt(mean((pred_c1 - log_lambda)^2)),
            rmse_ctrl2 = sqrt(mean((pred_c2 - log_lambda)^2)))

# write.csv(crossval_summary, 'analysis/stan_movwin_crossval_summary.csv', row.names = F)




##### cross validation via library 'caret' (in ML framework), for control models  ---------------------------------------------------------

CrossValML <- function(species_focal) {
  # getlambdas for focal focal species
  xx <- dplyr::filter(xx_full, species == unique(species_focal))
  
  # calculate monthly precipitation values
  years   <- xx$year %>% unique() %>% sort()
  precip_l  <- lapply(years, precip_form, precip_long, "precip")
  
  # put all in matrix form 
  x_precip  <- t(mat_form(precip_l, years))
  x_precip_means <- rowMeans(x_precip)   # averages over entire window (for control model #2)
  
  # add mean precip data (average over all 24 months) to xx
  xx$x_precip_means <- x_precip_means
  
  # add intercept column for caret::train() (newest version of caret gives warnings for intercept-only model)
  xx$intercept <- 1
  
  # perform loo cross-val via caret function train, for the ctrl1 and ctrl2 models
  mod_ctrl1 <- train(log_lambda ~ intercept, method = 'lm', data = xx, trControl = trainControl(method = 'LOOCV'))
  mod_ctrl2 <- train(log_lambda ~ x_precip_means, method = 'lm', data = xx, trControl = trainControl(method = 'LOOCV'))
  
  # extract out-of-sample rmse
  rmse_ctrl1_ml <- as.numeric(mod_ctrl1$results['RMSE'])
  rmse_ctrl2_ml <- as.numeric(mod_ctrl2$results['RMSE'])
  
  return(data.frame(rmse_ctrl1_ml, rmse_ctrl2_ml))
}

# perform leave-one-out cross validation for all 5 species, both control models
crossval_summary_ml <- xx_full %>%
  group_by(species) %>%
  do(CrossValML(species_focal = .$species)) %>% ungroup()

# write.csv(crossval_summary_ml, 'analysis/stan_movwin_crossval_summary_maxlik.csv', row.names = F)


# compare rmse from bayesian and ml cross-validation (for ctrl1 and ctrl2 only)
# values very similar, wohoo!
crossval_summary
crossval_summary_ml

