# first try to automate graphs 
setwd("C:/cloud/MEGA/projects/sAPROPOS/")
source("~/moving_windows/format_data.R")
library(dplyr)
library(tidyr)
library(testthat)
options(stringsAsFactors = F )


# read lambda/clim data ---------------------------------------------------------------------------------
lam     <- read.csv("C:/cloud/Dropbox/sAPROPOS project/DemogData/lambdas_6tr.csv", stringsAsFactors = F) 
m_info  <- read.csv("C:/cloud/MEGA/Projects/sApropos/MatrixEndMonth_information.csv", stringsAsFactors = F)
clim    <- read.csv("C:/cloud/Dropbox/sAPROPOS project/DemogData/precip_fc_demo.csv",  stringsAsFactors = F) 

# info on summaries ---------------------------------------------------------------------
res_folder<- "supercomp/res_7.26" 
sum_files <- list.files(res_folder)[grep("mod_summaries_", list.files(res_folder) )]


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



# single species graphs --------------------------------------------------------------------
post_files  <- list.files(res_folder)[grep("posterior_", list.files(res_folder) )]

# create species list based on what's available!
spp_list    <- intersect(gsub(".csv", "", gsub("mod_summaries_", "", sum_files) ),
                         gsub(".csv", "", gsub("posterior_", "", post_files) ) )
  

# function: plot gaussian models 
spp_plot_gaus <- function(spp_name){

  # format by species 
  m_back        <- 24     # months back
  
  # lambda data
  spp_lambdas   <- format_species(spp_name, lambdas)
  
  # climate data
  clim_separate <- clim_list(spp_name, clim, spp_lambdas)
  clim_detrnded <- lapply(clim_separate, clim_detrend)
  clim_mats     <- Map(clim_long, clim_detrnded, spp_lambdas, m_back)
  
  # model data
  mod_data          <- lambda_plus_clim(spp_lambdas, clim_mats)
  mod_data$climate  <- mod_data$climate #/ diff(range(mod_data$climate))
  
  
  # read mean values, whole posterior, GAUS MODELS ONLY
  post_df <- read.csv(paste0("C:/cloud/MEGA/Projects/sApropos/",res_folder,"/posterior_",
                      spp_name,".csv") ) %>%
                dplyr::select(-grep("log_lik", names(.)) ) %>%
                subset( model == "gaus" )
  mean_df <- read.csv(paste0("C:/cloud/MEGA/Projects/sApropos/",res_folder,"/mod_summaries_",
                      spp_name,".csv") ) %>%
                subset( model == "gaus" )
  
  
  # plots
  
  # posterior values and means
  pars_id   <- seq(1, (nrow(post_df) / 2), length.out = 100) %>% round
  par_smpl  <- post_df[pars_id,]
  
  # 
  month_x   <- seq(1,m_back,1)
  
  # 100 posterior samples of sensitivity
  store_norm <- function(i){
    sens <- dnorm(seq(1,m_back,1), par_smpl[i,"sens_mu"], par_smpl[i,"sens_sd"])
    sens / sum(sens)
  }
  # store 100 posterior samples of month sensitivity
  month_sens  <- do.call(rbind, lapply(c(1:100), store_norm) )
  
  # mean sensitivities sensitivities
  mean_m_sen  <- dnorm(seq(1, m_back,1), as.numeric(mean_df["sens_mu_mean"]), 
                       as.numeric(mean_df["sens_sd_mean"]) )
  mean_m_sen  <- mean_m_sen / sum(mean_m_sen)
  
  
  
  # figure: i) sensitivity to monthly precip; ii) log_lambda~climate ---------------------------
  tiff(paste0("results/plots/",spp_name,".tiff"),
       unit="in", width=6.3,height=3.5,res=600,compression="lzw")
  
  par(mfrow=c(1,2), mar=c(3.2,3.2,0.3,0.1), mgp = c(2,0.7,0), cex.lab = 1.2)
  
  # monthly sensitivities 
  plot(month_x, month_x,
       type="n", ylim = c(0,max(month_sens)),
       xlab = "Month", ylab = "Sensitivity to monthly precip.")
  
  # plot 100 posterior samples
  lapply(c(1:100), function(i) 
    lines(month_x, 
          month_sens[i,], 
          lwd=0.5, col = "grey") )
  # mean
  lines(month_x, mean_m_sen, lwd = 2)
  
  
  # prediction
  
  # climate predictor
  clim_mat   <- sweep(mod_data$climate, 2, mean_m_sen, "*") %>%
                  apply(1,sum)
  
  # plotting
  clim_x      <- stack(clim_mat) %>%
                    dplyr::select(-2) %>%
                    setNames(c("x_clim")) %>%
                    bind_cols(mod_data$lambdas)
  
  # plot 
  plot(log_lambda ~ x_clim, data = clim_x, pch = 16,
         xlab = "Precip. (weighted avg.)",
         ylab = expression(lambda))
  # plot posteriors
  lapply(1:100, function(i) 
    abline(a = par_smpl[i,"alpha"], col = "grey",
           b = par_smpl[i,"beta"], lwd = 2)
  )
  
  # plot means
  points(log_lambda ~ x_clim, data = clim_x, pch = 16,
       xlab = "Precip. (weighted avg.)",
       ylab = expression(lambda), new=T)
  abline(a = as.numeric(mean_df["alpha_mean"]),
         b = as.numeric(mean_df["beta_mean"]), lwd = 2)
  
  dev.off()

}

# produce the plots for 
lapply(spp_list, spp_plot_gaus)
