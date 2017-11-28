# graphs of effect sizes
setwd("C:/cloud/MEGA/projects/sAPROPOS/")
source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(testthat)
options(stringsAsFactors = F )

m_back        <- 24
clim_var      <- c("airt","precip")
spp_cat       <- read.csv("sApropos_spp_raunkiaer.csv") %>% rename( species = SpeciesAuthor )
sens_df       <- read.csv("allCOMPADRE_COMADREOutput.csv") %>%
                    #subset( MatrixComposite == "Individual" &
                    subset( MatrixComposite == "Mean" &
                            Treatment == "Unmanipulated" ) %>%
                    rename( species = SpeciesAuthor,
                            surv    = Ssurv,
                            grow    = Sgrow,
                            fec     = Srep) %>%
                    dplyr::select( species, surv, grow, fec) %>%
                    group_by( species ) %>%
                    summarize_all( mean, na.rm=T ) %>%
                    subset( !(is.nan(surv) & is.nan(grow) & is.nan(fec)) ) %>%
                    gather(response, sensitivity, -species )

# create the path of each summary file
path_summaries <- function(x,y){
  post_files <- Filter(function(x) grepl("mod_summaries",x), x)
  paste0(y,post_files)
}


# grab posterior mean directly
grab_summ_means <- function(x, response){
  
  # extract species name
  spp_n     <- gsub( paste0("results/moving_windows/",response,"/[a-z]{4,6}/mod_summaries_"),"",x)
  spp_name  <- gsub(".csv","",spp_n)
  
  # extract climate variable name
  # substitute 
  # 1 everything before "/posterior_"
  # 2 everything after  "response/"
  clim_var  <- gsub(paste0("^(.*?)/",response,"/|/mod_summaries_?.*"),"",x)
  
  mean_by_spp <- read.csv( x ) %>%
    dplyr::select( model,alpha_mean,beta_mean ) %>%
    as.data.frame %>%
    mutate(species  = spp_name,
           clim_var = clim_var )
  
}

# get means from summary files
extract_summ_means <- function(response, clim_var, mod_name = "ctrl2" ){
  
  # set up list of paths go get files
  wds           <- paste0("results/moving_windows/",response,"/",clim_var,"/") %>% as.list
  raw_file_vec  <- lapply(wds, list.files)
  # paths of files containing summaries
  summ_paths_l  <- Map(path_summaries, raw_file_vec, wds)
  summ_paths    <- Reduce(function(...) c(...), summ_paths_l) %>% as.list
  
  # actually get summary means
  means_l       <- lapply(summ_paths, grab_summ_means, response)
  means_out     <- Reduce(function(...) bind_rows(...), means_l) %>%
                      mutate( response = response) %>%
                      subset( model    == mod_name) %>%
                      rename_( .dots = setNames("beta_mean",  paste0("beta_",response)) ) %>%
                      rename_( .dots = setNames("alpha_mean", paste0("alpha_",response)) ) %>%
                      dplyr::select( -response )
  
  return(means_out)

}

# read and format beta values
means_s_pre <- extract_summ_means("SurvSSDPreRep",  c("airt","precip") )
means_s_rep <- extract_summ_means("SurvSSDRep",     c("airt","precip") )
means_g_pre <- extract_summ_means("GrowSSDPreRep",  c("airt","precip") )
means_g_rep <- extract_summ_means("GrowSSDRep",     c("airt","precip") )

# substitute string with another string
subs_with   <- function(x, orig, new){
  
  replace(x, x == orig, new)
  
}

# put betas all in one data frame
means_beta   <- Reduce( function(...) full_join(...), 
                       list(means_s_pre, means_s_rep,
                            means_g_pre, means_g_rep) ) %>%
                  dplyr::select(species, clim_var, 
                                beta_SurvSSDPreRep, beta_SurvSSDRep, 
                                beta_GrowSSDPreRep, beta_GrowSSDRep ) %>%
                  gather( response, beta, -species, -clim_var ) %>%
                  mutate( response = subs_with(response, "beta_SurvSSDPreRep", "survPreRep") ) %>%
                  mutate( response = subs_with(response, "beta_SurvSSDRep", "survRep") ) %>%
                  mutate( response = subs_with(response, "beta_GrowSSDPreRep", "growPreRep") ) %>%
                  mutate( response = subs_with(response, "beta_GrowSSDRep", "growRep") ) 

# put betas all in one data frame
means_alph    <- Reduce( function(...) full_join(...), 
                       list(means_s_pre, means_s_rep,
                            means_g_pre, means_g_rep) ) %>%
                  dplyr::select(species, clim_var, 
                                alpha_SurvSSDPreRep, alpha_SurvSSDRep, 
                                alpha_GrowSSDPreRep, alpha_GrowSSDRep ) %>%
                  gather( response, alpha, -species, -clim_var ) %>%
                  mutate( response = subs_with(response, "alpha_SurvSSDPreRep", "survPreRep") ) %>%
                  mutate( response = subs_with(response, "alpha_SurvSSDRep", "survRep") ) %>%
                  mutate( response = subs_with(response, "alpha_GrowSSDPreRep", "growPreRep") ) %>%
                  mutate( response = subs_with(response, "alpha_GrowSSDRep", "growRep") ) 

# means together
means_all    <- full_join(means_alph, means_beta)

# sensitivities 
sens_ssd       <- read.csv("allCOMPADRE_COMADREOutput.csv") %>%
                    dplyr::select(SpeciesAuthor,
                                  SsurvPreRep, SsurvRep,
                                  SgrowPreRep, SgrowRep ) %>%
                    rename( species = SpeciesAuthor) %>%
                    group_by( species ) %>%
                    summarize_all( mean, na.rm=T ) %>%
                    subset( !(is.nan(SsurvPreRep) & is.nan(SsurvRep) & 
                              is.nan(SgrowPreRep) & is.nan(SgrowRep)) ) %>%
                    gather(response, sensitivity, -species ) %>%
                    mutate( response = subs_with(response, "SsurvPreRep", "survPreRep") ) %>%
                    mutate( response = subs_with(response, "SsurvRep", "survRep") ) %>%
                    mutate( response = subs_with(response, "SgrowPreRep", "growPreRep") ) %>%
                    mutate( response = subs_with(response, "SgrowRep", "growRep") ) 

# standardized betas
means_sens     <- inner_join(means_all, sens_ssd) %>%
                    mutate( clim_sens        = ( beta * exp(alpha) ) / ( (exp(alpha) +1)^2 ) ) %>%
                    mutate( clim_sens        = abs(clim_sens) ) %>%
                    mutate( clim_sens_scaled = sensitivity * clim_sens) %>%
                    dplyr::select(-alpha,-beta,-sensitivity,-clim_sens) %>%
                    spread(species, clim_sens_scaled)


# plot -----------------------------------------------------------------------------------------------

# separate by weather predictor
mean_airt      <- subset(means_sens, clim_var == "airt") %>% .[c(3,1,4,2),]
mean_prec      <- subset(means_sens, clim_var == "precip") %>% .[c(3,1,4,2),]

# set up color coding 
col_df <- data.frame( species = names(mean_airt)[-c(1,2)]) %>%
            left_join( spp_cat ) %>%
            mutate( col = "blue") %>%
            mutate( col = replace( col,  
                                   GrowthFormRaunkiaer == "Geophyte" | 
                                   GrowthFormRaunkiaer == "Chamaephyte",
                                   "brown" ) )


# sensitivity by vr and reproductive stage
tiff("results/moving_windows/sens_by_repr_vr.tiff",
     unit="in", width=4.5, height=6.3, res=400, compression="lzw")

par( mfrow = c(2,1), mar = c(3,4,0.1,1.5), mgp = c(2,0.7,0) )

matplot( mean_airt[,-c(1,2)], cex.lab=1.1, lwd = 2, ylim = c(0,1),
         col = col_df$col, lty = 1, type = "l", xaxt = "n", 
         ylab = expression("Sensitivity of  "*lambda*"  to temperature") )
axis(1, at = c(1,2,3,4),
     labels = c("Surv_Pre","Grow_Pre","Surv_Rep","Grow_Rep") )

matplot( mean_prec[,-c(1,2)], cex.lab=1.1, lwd = 2, ylim = c(0,1),
         col = col_df$col, lty = 1, type = "l", xaxt = "n",
         ylab = expression("Sensitivity of "*lambda*" to precipitation") )
axis(1, at = c(1,2,3,4), 
     labels = c("Surv_Pre","Grow_Pre","Surv_Rep","Grow_Rep") )

dev.off()
