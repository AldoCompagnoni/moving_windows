#bjtwd<-"C:/Users/admin_bjt162/Dropbox/A.Current/Ongoing_Collab_Research/sApropos project/"
rm(list=ls())
setwd("C:/cloud/MEGA/Projects/sApropos/")
source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(dismo)
library(mgcv)
library(testthat)

# climate predictor, months back, max. number of knots
clim_var  <- "airt"
m_back    <- 24    
knots     <- 3


# read data -----------------------------------------------------------------------------------------
lam       <- read.csv("lambdas_6tr.csv", stringsAsFactors = F) 
m_info    <- read.csv("MatrixEndMonth_information.csv", stringsAsFactors = F)
clim_fc   <- data.table::fread(paste0(clim_var,"_fc_demo.csv"),  stringsAsFactors = F)
clim_35   <- read.csv( paste0("monthly_",clim_var,"_Dalgleish.csv") )
clim      <- list(clim_fc, clim_35)
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
  unique_yr         <- mod_data$lambdas$year %>% unique
  
  # format data for splines --------------------------------------------------------
  
  # precipitation matrix
  pmat      <- mod_data$climate %>% setNames(NULL) %>% as.matrix
  pmean     <- apply(pmat, 1, mean, na.rm = T)
  
  # set up lags
  lags      <- matrix(0, nrow(mod_data$climate), ncol(mod_data$climate))
  for(i in 1:ncol(lags)) lags[,i]=i
  lagsm     <- as.matrix(lags)
  
  # assemble final data frame
  # pmat <- pmat/diff(range(pmat));
  dat       <- cbind(mod_data$lambdas, pmat)
  dat       <- as.data.frame(dat, colnames=TRUE)
  dat$lags  <- lags
  dat$pmat  <- pmat
  dat$pmean <- pmean
  
  # model fits -----------------------------------------------------------------------------------
  
  # fit full model
  if( length(unique(dat$population)) > 1 ){
    mod_full  <-gam(log_lambda ~ s(lags, k=knots, by=pmat, bs="cs"), # + population,
                    data=dat,method="GCV.Cp",gamma=1.4, na.action = na.omit)
  }else{
    mod_full  <-gam(log_lambda ~ s(lags, k=knots, by=pmat, bs="cs"),
                    data=dat,method="GCV.Cp",gamma=1.4, na.action = na.omit)
  }
  
  # crossvalidation function
  crxval_spline <- function(i, dat, pmat){
    
    # select excluded year
    yr          <- unique_yr[i]
    
    # train and test sets, null prediction
    train_set   <- subset( dat, year != yr )
    test_set    <- subset( dat, year == yr )
    
    # model: leave-one-year-out
    if( length(unique(dat$population)) > 1 ){
      mod_spline  <- gam(log_lambda ~ s(lags, k=knots, by=pmat, bs="cs",sp=mod_full$sp),# + population,
                         method="GCV.Cp",gamma=1.4, na.action = na.omit, data=train_set)
      mod_lm      <- lm(log_lambda ~ pmean, data=train_set)# + population
    }else{
      mod_spline  <- gam(log_lambda ~ s(lags, k=knots, by=pmat, bs="cs",sp=mod_full$sp),
                         method="GCV.Cp",gamma=1.4, na.action = na.omit, data=train_set)
      mod_lm      <- lm(log_lambda ~ pmean, data=train_set)
    }
    
    # predictions leave-one-out: spline, linear model
    pred_null   <- mean( train_set$log_lambda )
    pred_lm     <- predict(mod_lm,     newdata = test_set, type="response")
    pred_spline <- predict(mod_spline, newdata = test_set, type="response")
    
    # data frame of predictions and "design" (year-by-population)
    pred_df     <- data.frame(pred_null   = pred_null, 
                              pred_lm     = pred_lm,
                              pred_spline = pred_spline,
                              stringsAsFactors = F
                              )
    design_pred <- dplyr::select(test_set, year, population)
    
    bind_cols(design_pred, pred_df)
    
  }
  
  # apply function across samples, and create data frame
  crxval_l    <- lapply(1:length(unique_yr), crxval_spline, dat, pmat)
  crxval_df   <- Reduce(function(...) rbind(...), crxval_l) %>%
                    left_join( dplyr::select(dat,year, population,log_lambda) ) # %>%
                    # subset( !is.na(pred_spline) )
  
  # Mean squared error
  dev0        <- calc.deviance(crxval_df$log_lambda, crxval_df$pred_null, 
                               weights = rep(1, nrow(crxval_df) ),  
                               family="gaussian", calc.mean = TRUE)  
  dev1        <- calc.deviance(crxval_df$log_lambda, crxval_df$pred_lm, 
                               weights = rep(1, nrow(crxval_df) ),  
                               family="gaussian", calc.mean = TRUE)
  dev2        <- calc.deviance(crxval_df$log_lambda, crxval_df$pred_spline,
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
  pops <- crxval_df$population %>% unique
  
  # all lambdas (for limits of y-axis)
  all_llam  <- c(crxval_df$log_lambda, 
                 crxval_df$pred_null, 
                 crxval_df$pred_lm, 
                 crxval_df$pred_spline) 
  
  # graph lambdas 
  for(p in 1:length(pops)){
    
    # population index
    pop_i       <- which(crxval_df$population == pops[p])
    
    # prepare data for plotting
    plot_data   <- crxval_df[pop_i,]
    plot_null   <- plot_data$pred_null[pop_i]
    plot_lm     <- plot_data$pred_lm[pop_i]
    plot_spline <- plot_data$pred_spline[pop_i]
    
    if( p == 1){
      plot(log_lambda ~ year, type="o", ylim = range(all_llam),
           ylab="logLambda", xlab="Left-out year", col="black", data = plot_data)
      
    } else{
      points(plot_data$year, plot_data$log_lambda, type="o",
           ylab="logLambda", xlab="Left-out year", col="black")
    }
    points(plot_data$year, plot_data$pred_null, type="o", col="red", lty = p)
    points(plot_data$year, plot_data$pred_lm, type="o", col="brown", lty = p)
    points(plot_data$year, plot_data$pred_spline, type="o", col="blue", lty = p)
    
  }
  
  legend("topright", c("Data", "NULL", "LM", "Spline"),
         col=c("black", "red", "brown", "blue"),
         lty = c(1, 1, 1), cex = 0.8, lwd = 2 , bty="n")
  
  # legends deviance
  dev.0=round(dev0,5)
  dev.1=round(dev1,5)
  dev.2=round(dev2,5)
  
  legend("bottomleft", 
         c(paste0("NULL mse: ",dev.0), 
           paste0("LM           : ",dev.1), 
           paste0("Spline       : ",dev.2)),
         col=c( rep("black",length(pops)+1), "red", "blue"),
         lty = c(seq_along(pops), 1, 1, 1), bty="n", cex = 0.8)
  
  dev.off()

  mod_sum_l[[ii]] <- crxval_df %>%
                        mutate( dev0 = dev0,
                                dev1 = dev1,
                                dev2 = dev2,
                                species = spp_name )
                    
}

# summary figures
mod_sum_df <- Reduce(function(...) rbind(...), mod_sum_l)
write.csv(mod_sum_df, 
          paste0("results/splines/loyo/summaries/spline_",
                 clim_var,m_back,"_summaries.csv"),
          row.names=F)
