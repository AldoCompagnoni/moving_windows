#bjtwd<-"C:/Users/admin_bjt162/Dropbox/A.Current/Ongoing_Collab_Research/sApropos project/"
setwd("C:/cloud/MEGA/Projects/sApropos/")
source("~/moving_windows/format_data.R")
library(tidyverse)
library(dismo)
library(mgcv)
library(testhtat)

# climate predictor, months back, max. number of knots
clim_var  <- "pet"
m_back    <- 24    
knots     <- 8

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv("lambdas_6tr.csv", stringsAsFactors = F) 
m_info    <- read.csv("MatrixEndMonth_information.csv", stringsAsFactors = F)
clim      <- read.csv(paste0(clim_var,"_fc_demo.csv"),  stringsAsFactors = F)
spp       <- clim$species %>% unique


# add monthg info to lambda information
month_add <- m_info %>%
                mutate(SpeciesAuthor = trimws(SpeciesAuthor) ) %>%
                dplyr::select(SpeciesAuthor, MatrixEndMonth)
lam_add   <- subset(lam, SpeciesAuthor %in% month_add$SpeciesAuthor) %>%  
                dplyr::select(-MatrixEndMonth) %>%
                inner_join(month_add)
lam_min   <- subset(lam, SpeciesAuthor %in% setdiff(lam$SpeciesAuthor, m_info$SpeciesAuthor) )
lambdas   <- bind_rows(lam_min, lam_add) %>%
                subset( !is.na(MatrixEndMonth) )

mod_sum_l <- list()
spp_list  <- lambdas$SpeciesAuthor %>% unique   
spp_list  <- spp_list[ -c(2,3,4,8,9,16,21,22,24) ] #,3,8,14,18,19,20)]

# run models and store pictures ------------------------------------------------------
for(ii in 1:length(spp_list)){
  
  # set species
  spp_name      <- spp_list[ii] 
  # lambda data
  spp_lambdas   <- format_species(spp_name, lambdas)
  
  # climate data
  clim_separate <- clim_list(spp_name, clim, spp_lambdas)
  clim_detrnded <- lapply(clim_separate, clim_detrend, clim_var)
  clim_mats     <- Map(clim_long, clim_detrnded, spp_lambdas, m_back)
  
  # model data
  mod_data          <- lambda_plus_clim(spp_lambdas, clim_mats)
  mod_data$climate  <- mod_data$climate #/ diff(range(mod_data$climate))
  
  # tests
  expect_equal(length(spp_lambdas), length(clim_separate) )
  
  # unique years, the basis of crossvalidation samples
  unique_yr         <- arrange(mod_data$lambdas,year) %>% .$year %>% unique
  
  # format data for splines --------------------------------------------------------
  
  # precipitation matrix
  pmat      <- mod_data$climate %>% setNames(NULL) %>% as.matrix
  
  # set up lags
  lags      <- matrix(0, nrow(mod_data$climate), ncol(mod_data$climate)); 
  for(i in 1:ncol(lags)) lags[,i]=i; 
  lagsm     <- as.matrix(lags); 
  
  # assemble final data frame
  # pmat <- pmat/diff(range(pmat));
  dat       <- cbind(mod_data$lambdas, pmat)
  dat       <- as.data.frame(dat, colnames=TRUE)
  dat$lags  <- lags
  dat$pmat  <- pmat
  
  
  # model fits -----------------------------------------------------------------------------------
  
  # fit full model
  if( length(unique(dat$population)) > 1 ){
    mod_full  <-gam(log_lambda ~ s(lags, k=knots, by=pmat, bs="cs"), # + population,
                    data=dat,method="GCV.Cp",gamma=1.4)
  } else{
    mod_full  <-gam(log_lambda ~ s(lags, k=knots, by=pmat, bs="cs"),
                    data=dat,method="GCV.Cp",gamma=1.4)
  } 
  
  # crossvalidation function
  crxval_spline <- function(i, dat, pmat){
    
    # select excluded year
    yr          <- unique_yr[i]
    
    # train and test sets, null prediction
    train_set   <- subset( dat, year != yr )
    test_set    <- subset( dat, year == yr )
    pred_null   <- mean( train_set$log_lambda )
    
    # model: leave-one-out
    if( length(unique(dat$population)) > 1 ){
      mod_loo     <- gam(log_lambda ~ s(lags, k=knots, by=pmat, bs="cs",sp=mod_full$sp),# + population,
                         method="GCV.Cp",gamma=1.4, data=train_set)
    }else{
      mod_loo     <- gam(log_lambda ~ s(lags, k=knots, by=pmat, bs="cs",sp=mod_full$sp),
                         method="GCV.Cp",gamma=1.4, data=train_set)
    }
    
    
    # prediction: leave-one-out
    pred_oo     <- predict(mod_loo, newdata = test_set, type="response")
    
    data.frame( pred_null = pred_null, 
                pred_oo   = pred_oo )
    
  }
  
  # apply function across samples, and create data frame
  crxval_l    <- lapply(1:length(unique_yr), crxval_spline, dat, pmat)
  crxval_df   <- Reduce(function(...) rbind(...), crxval_l)
  
  # Mean squared error
  dev0        <- calc.deviance(crxval_df$pred_null, dat$log_lambda, 
                               weights = rep(1, nrow(crxval_df) ),  
                               family="gaussian", calc.mean = TRUE)  
  dev1        <- calc.deviance(crxval_df$pred_oo, dat$log_lambda, 
                               weights = rep(1, nrow(crxval_df) ),  
                               family="gaussian", calc.mean = TRUE)
  
 
  # plot results ---------------------------------------------------------------
  tiff(paste0("results/splines/loyo/plots/",clim_var,m_back,"/",spp_name,".tiff"),
       unit="in", width=6.3, height=3.15, res=400, compression="lzw")
  
  par(mfrow=c(1,2), mar = c(3,3,1.7,0.1), mgp = c(1.8,0.7,0))
  plot(mod_full)
  mtext(spp_name, side = 3, line = 0.4,
        at = par("usr")[2], cex = 1.5)
  abline(h=0)
  
  # plot populations separately
  pops <- dat$population %>% unique
  
  # all lambdas (for limits of y-axis)
  all_llam  <- c(dat$log_lambda, crxval_df$pred_null, crxval_df$pred_oo) 
  
  # graph lambdas 
  for(p in 1:length(pops)){
    
    # population index
    pop_i     <- which(dat$population == pops[p])
    
    # prepare data for plotting
    plot_data <- dat[pop_i,]
    plot_null <- crxval_df$pred_null[pop_i]
    plot_pred <- crxval_df$pred_oo[pop_i]
    
    if( p == 1){
      plot(plot_data$log_lambda ~ plot_data$year, type="o", ylim = range(all_llam),
           ylab="logLambda", xlab="Left-out year", col="black")
      
    } else{
      points(plot_data$log_lambda ~ plot_data$year, type="o",
           ylab="logLambda", xlab="Left-out year", col="black")
    }
    points(plot_null ~ plot_data$year, type="o", col="red", lty = p)
    points(plot_pred ~ plot_data$year, type="o", col="blue", lty = p)
    
  }
  
  legend("topright", c("Data", "NULL", "crossval"),
         col=c("black", "red", "blue"),
         lty = c(1, 1, 1), cex = 0.8, lwd = 2 , bty="n")
  
  # legends deviance
  dev.0=round(dev0,5)
  dev.1=round(dev1,5)
  
  legend("bottomleft", 
         c(paste0("NULL mse: ",dev.0), paste0("Full. Mod. mse: ",dev.1)), 
         col=c( rep("black",length(pops)+1), "red", "blue"),
         lty = c(seq_along(pops), 1, 1, 1), bty="n", cex = 0.8)
  
  dev.off()

  mod_sum_l[[ii]] <- crxval_df %>%
                        mutate( dev0 = dev0,
                                dev1 = dev1,
                                species = spp_name )
                    
}

# summary figures
mod_sum_df <- Reduce(function(...) rbind(...), mod_sum_l)
write.csv(mod_sum_df, 
          paste0("results/splines/loyo/summaries/spline_",
                 clim_var,m_back,"_summaries.csv"),
          row.names=F)
