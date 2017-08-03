#bjtwd<-"C:/Users/admin_bjt162/Dropbox/A.Current/Ongoing_Collab_Research/sApropos project/"
setwd("C:/cloud/Dropbox/sAPROPOS project/")
source("~/moving_windows/format_data.R")
library(tidyverse)
library(dismo)
library(mgcv)


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
lam_add   <- subset(lam, SpeciesAuthor %in% month_add$SpeciesAuthor) %>%  
                dplyr::select(-MatrixEndMonth) %>%
                inner_join(month_add)
lam_min   <- subset(lam, SpeciesAuthor %in% setdiff(lam$SpeciesAuthor, m_info$SpeciesAuthor) )
lambdas   <- bind_rows(lam_min, lam_add) %>%
                subset( !is.na(MatrixEndMonth) )

mod_sum_l <- list()
spp_list  <- lambdas$SpeciesAuthor %>% unique   
spp_list  <- spp_list[-c(2,3,4,8,9,16,21,22,24)] #,3,8,14,18,19,20)]

# perliminary format of species-level data ------------------------------------------------------
for(ii in 1:length(spp_list)){
  
  spp_name      <- spp_list[ii] # test run w/ spp number 1
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
  mod_full  <-gam(log_lambda ~ s(lags, k=10, by=pmat, bs="cs"),
                  data=dat,method="GCV.Cp",gamma=1.4) 
  
  # crossvalidation function
  crxval_spline <- function(i, dat, pmat){
    
    dati        <- dat[-i,]
    pred_null    <- mean( dat$log_lambda[-i] )
    
    mod_loo     <- gam(log_lambda ~ s(lags, k=10, by=pmat, bs="cs",sp=mod_full$sp),
                       method="GCV.Cp",gamma=1.4, data=dati) 
    
    pred_oo     <- predict(mod_loo, newdata = dat[i,], type="response")
    
    data.frame( pred_null = pred_null, 
                pred_oo   = pred_oo )
    
  }
  # apply function, and create data frame
  crxval_l    <- lapply(1:nrow(dat), crxval_spline, dat, pmat)
  crxval_df   <- Reduce(function(...) rbind(...), crxval_l)
  
  # Mean squared error
  dev0        <- calc.deviance(crxval_df$pred_null, dat$log_lambda, 
                               weights = rep(1, nrow(crxval_df) ),  
                               family="gaussian", calc.mean = TRUE)  
  dev1        <- calc.deviance(crxval_df$pred_oo, dat$log_lambda, 
                               weights = rep(1, nrow(crxval_df) ),  
                               family="gaussian", calc.mean = TRUE)
  
 
  # plot results ---------------------------------------------------------------
  tiff(paste0("C:/cloud/MEGA/Projects/sApropos/results/splines/plots/",spp_name,".tiff"),
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
  
  legend("topright", c("full model", "NULL", "crossval"),
         col=c("black", "red", "blue"),
         lty = c(1, 1, 1), cex = 0.8, lwd = 2 , bty="n")
  
  # legends deviance
  dev.0=round(dev0,3)
  dev.1=round(dev1,3)
  
  legend("bottomleft", 
         c(paste0("NULL mse: ",dev.0), paste0("Full. Mod. mse: ",dev.1)), 
         col=c( rep("black",length(pops)+1), "red", "blue"),
         lty = c(seq_along(pops), 1, 1, 1), bty="n", cex = 0.8)
  
  dev.off()

  mod_sum_l[[ii]] <- crxval_df %>%
                        mutate( dev0 = dev0,
                                dev1 = dev1,
                                species = spp_name)
                    
}

# summary figures
mod_sum_df <- Reduce(function(...) rbind(...), mod_sum_l) 
write.csv(mod_sum_df, 
          "C:/cloud/MEGA/Projects/sApropos/results/splines/spline_summaries.csv",
          row.names=F)
