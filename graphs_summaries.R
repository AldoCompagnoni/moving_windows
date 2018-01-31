# produce summary graphs
setwd("C:/cloud/MEGA/projects/sAPROPOS/")
source("C:/CODE/moving_windows/format_data_5.R")
library(tidyverse)
library(magrittr)
library(testthat)
options(stringsAsFactors = F )

m_back        <- 36
interval      <- NULL

# pipe-able Reduce_rbind function 
rbind_l     <- function(df_list){ Reduce(function(...) rbind(...), df_list) }

# summarize moving windows results by climate variable
summ_by_clim <- function(ii){
  
  clim_var <- input_df$clim_var[ii]
  response <- input_df$response[ii]
  interval <- input_df$interval[ii]
  
  # read lambda/clim data ---------------------------------------------------------------------------------
  lam     <- read.csv("all_demog_6tr.csv", stringsAsFactors = F) %>%
                subset( SpeciesAuthor != "Purshia_subintegra" )
  m_info  <- read.csv("MatrixEndMonth_information.csv", stringsAsFactors = F)
  clim    <- data.table::fread(paste0(clim_var,"_fc_hays.csv"),  stringsAsFactors = F)
  
  # add monthg info to lambda information
  month_add <- m_info %>%
                  mutate(SpeciesAuthor = trimws(SpeciesAuthor) ) %>%
                  dplyr::select(SpeciesAuthor, MatrixEndMonth)
  lam_add   <- subset(lam, SpeciesAuthor %in% month_add$SpeciesAuthor) %>%  
                  dplyr::select(-MatrixEndMonth) %>%
                  inner_join(month_add)
  lam_min   <- subset(lam, SpeciesAuthor %in% setdiff(lam$SpeciesAuthor, m_info$SpeciesAuthor) )
  lambdas   <- bind_rows(lam_min, lam_add) %>%
                  subset( !is.na(MatrixEndMonth) ) %>%
                  subset( SpeciesAuthor != "Eryngium_alpinum" )
  
  # summaries organism type/biome
  by_spp_rep <- function(var){
    
    lambdas %>%
      dplyr::select( c("SpeciesAuthor", var) ) %>%
      distinct %>%
      group_by_(var) %>%
      summarise(spp_n = n())
    
  }
  
  # replication of categories
  categ_rep <- lapply(c("OrganismType","Ecoregion","Order", "Class", "DicotMonoc", "Altitude"), 
                      by_spp_rep)
  
  # I chose to use the following categories
  categ     <- lambdas %>%
                  dplyr::select( c("SpeciesAuthor","OrganismType","Ecoregion", "DicotMonoc", "Class") ) %>%
                  rename( species = SpeciesAuthor) %>%
                  unique
  
  # summary info ----------------------------------------------------------------------------
  res_folder<- paste0("results/moving_windows/",response,"/",clim_var,interval) 
  sum_files <- list.files(res_folder)[grep("mod_summaries_", list.files(res_folder) )]
  crx_files <- list.files(res_folder)[grep("crossval_", list.files(res_folder) )]
  mod_summ  <- lapply(sum_files, function(x) read.csv(paste0(res_folder,"/",x)) ) %>%
                  setNames( gsub("mod_summaries_", "", sum_files ) ) %>%
                  setNames( gsub(".csv", "", names(.) ) )
  replic    <- lapply(crx_files, function(x) read.csv(paste0(res_folder,"/",x)) %>% nrow ) %>%
                  setNames( gsub("crossval_pred_diag_", "", crx_files ) ) %>%
                  setNames( gsub(".csv", "", names(.) ) ) %>%
                  unlist %>% t %>% t %>% as.data.frame %>% 
                  tibble::rownames_to_column(var = "species") %>%
                  rename( rep_n = V1)
  if( clim_var != "soilmoist"){
    mod_splin <- read.csv( paste0("results/splines/summaries/spline_",clim_var,"24_summaries.csv") ) %>%
                    dplyr::select(species,dev0,dev1) %>%
                    unique %>%
                    mutate( model_climate = 0 ) %>%
                    mutate( model_climate = replace( model_climate, which(dev0 > dev1), 1) ) %>%
                    setNames( c("species", "dev0_spline", "dev1spline", "model_climate_spline") )
  }
  
  # predictive accuracy by measure of fit
  pred_acc_by_mof <- function(mof){
    
    if( mof == "mse"){
      tmp <- lapply(mod_summ, function(x) arrange_(x, mof)) %>%
                lapply(function(x) x[1,c(mof,"model")]) 
    }
    if( mof == "elpd"){
      tmp <- lapply(mod_summ, function(x) x[order(x[,mof], decreasing=T),] ) %>%
                lapply(function(x) x[1,c(mof,"model")]) 
    }
    tmp <- Map(function(x,y) tibble::add_column(x, species = y, .before = 1), 
               tmp, names(tmp) )
    
    Reduce(function(...) rbind(...), tmp) %>%
      inner_join( replic ) %>%
      setNames( c("species", mof , paste0("model_", mof), paste0("rep_n_", mof)) )
    
  }
  pred_acc_l  <- Map(pred_acc_by_mof, c("mse", "elpd") )
  
  # replication in years
  rep_yr <- function(spp_name){
    
    spp_df    <- subset(lambdas, SpeciesAuthor == spp_name)
    spp_yr    <- spp_df$MatrixEndYear %>% unique %>% length
    
    data.frame( species = spp_name,
                rep_yr  = spp_yr)
    
  }
  spp_rep_yr  <- lapply(unique(lambdas$SpeciesAuthor), rep_yr) %>% rbind_l
  
  # convergence info
  converg_df  <- Map(function(x,y) dplyr::select(x,model,rhat_high) %>% mutate( species = y), 
                     mod_summ, names(mod_summ) ) %>% rbind_l
  
  # set up categories for summary graphs
  mod_climate <- Reduce(function(...) merge(...), pred_acc_l) %>%
                    rename( rep_n = rep_n_mse) %>%
                    left_join( spp_rep_yr ) %>%
                    left_join(converg_df) %>%
                    inner_join( categ ) %>%
                    mutate( model_climate = sub("ctrl1", "0", model_mse)) %>%
                    mutate( model_climate = sub("ctrl2|expp|gaus", "1", model_climate)) %>%
                    mutate( model_climate = as.numeric(model_climate)) %>%
                    mutate( Ecoregion = as.factor(Ecoregion)) %>%
                    mutate( DicotMonoc = as.factor(DicotMonoc)) %>%
                    mutate( Class = as.factor(Class) ) %>%
                    mutate( OrganismType = as.factor(OrganismType) ) %>%
                    mutate( response_var = response,
                            clim_variable = clim_var )
  
  return(mod_climate)
  
}
# # put it all together
# mw_summ_l   <- lapply(c("precip", "pet", "airt", "airt_gdd"), summ_by_clim )
# # Also: color code the data frame 
# mw_summ_df  <- Reduce(function(...) rbind(...), mw_summ_l) %>%
#                   merge( data.frame( clim_var = c("precip", "pet", "airt", "airt_gdd"),
#                                      color    = c("blue", "orange", "red", "green") ) ) %>%
#                   mutate(method = "mov_win")
# put it all together
# mw_summ_l   <- lapply(c("precip", "airt", "soilmoist"), summ_by_clim, response)
# mw_summ_df  <- Reduce(function(...) rbind(...), mw_summ_l) %>%
#                   merge( data.frame( clim_var = c("precip", "airt",  "soilmoist"),
#                                      color    = c("blue", "red", "brown") ) ) %>%
#                   mutate(method = "mov_win")
input_df    <- expand.grid( clim_var = c("precip","airt"),
                            response = c("surv","grow","fec"),
                            interval = "" )

# ALL model information
mw_summ_df  <- lapply(1:nrow(input_df), summ_by_clim) %>%
                  rbind_l %>%
                  mutate(method = "mov_win")

# # convergence issues
# mw_summ_df %>% 
#   subset(model=="gaus") %>% 
#   subset( rhat_high > 0) %>% 
#   group_by(clim_variable, response_var) %>%
#   summarise( convergence_issues = n())
                  
                  
# add splines to the mix
read_spline_summ <- function(clim_var, m_back){
  return(
    read.csv( paste0("results/splines/summaries/spline_",clim_var,m_back,"_summaries.csv") ) %>%
      dplyr::select(species,dev0,dev1,dev2) %>%
      unique %>%
      setNames( c("species", "dev_NULL", "dev_lm",  "dev_spline") ) %>%
      mutate( model_climate = 0 ) %>%
      mutate( model_climate = replace( model_climate, which(dev_NULL > dev_lm | dev_NULL > dev_spline ), 1) ) %>%
      mutate( model_mse = "NULL") %>%
      mutate( model_mse = replace( model_mse, which(model_climate == 1 ), "spline") ) %>%
      mutate( model_mse = replace( model_mse, which(model_mse ==  "spline" & dev_spline > dev_lm ), "24mon") ) %>%
      mutate( clim_var = clim_var )
  )
}
spline_summ_l   <- lapply( c("airt","airt_gdd","pet","precip"), read_spline_summ, 12)
spline_summ_tmp <- Reduce( function(...) rbind(...), spline_summ_l ) %>%
                      merge( data.frame( clim_var = c("precip", "pet", "airt",  "airt_gdd"),
                                         color    = c("blue", "orange", "red", "green") ) ) %>%
                      mutate( method = "spline") %>%
                      rowwise %>%
                      mutate(mse = min(dev_NULL,dev_lm,dev_spline) )

# add information on replication, and climate sampled
rep_clim_info   <- dplyr::select(mw_summ_df, species, clim_var, rep_n, rep_yr,
                                 prop_rang, prop_yrs,prop_var,prop_var_r, mean_dev, mean_clim ) %>%
                        unique
spline_summ_df  <- left_join( spline_summ_tmp, rep_clim_info )


# summary tables ------------------------------------------------------------------------------------------------

# best models: AGGREGATE results
best_mod_agg <- dplyr::select(mw_summ_df, species, response_var, clim_variable, model_elpd) %>%
                  unique %>%
                  group_by( model_elpd ) %>%
                  summarise( best_mod_n = n() ) %>%
                  mutate( percent = round( best_mod_n/sum(best_mod_n), 2)*100 )
write.csv(best_mod_agg, "results/moving_windows/best_model_vital_rates.csv", row.names=F)

# set up barplot data frame
best_mod_all <- dplyr::select(mw_summ_df, species, response_var, clim_variable, model_elpd) %>%
                  unique %>%
                  group_by( response_var, clim_variable, model_elpd ) %>%
                  summarise( best_mod_n = n() ) %>%
                  mutate( percent = round( best_mod_n/sum(best_mod_n), 2)*100 ) %>%
                  arrange(model_elpd, response_var, clim_variable)
# all cases
expand_cases <- expand.grid( model_elpd    = c("ctrl1","ctrl2","expp","gaus", "gev", "simpl"),
                             response_var  = c("surv","grow","fec"),
                             clim_variable = c("airt","precip") )

# format model selection results to present in a table
form_mod_sel <- function( clim_var ){
  
  out <- expand_cases %>%
            subset( clim_variable == clim_var ) %>%
            left_join( subset(best_mod_all, clim_variable == clim_var) ) %>%
            mutate( best_mod_n = replace(best_mod_n, is.na(best_mod_n), 0),
                    percent    = replace(percent, is.na(percent), 0) ) %>%
            arrange(model_elpd, response_var)
  return(out)
  
}

# prepare model selection table
best_mod_all <- lapply(c("precip","airt"), form_mod_sel) %>% 
                  rbind_l %>%
                  mutate( response_var = as.character(response_var) ) %>%
                  mutate( response_var = replace(response_var, response_var == "surv", "Survival"),
                          response_var = replace(response_var, response_var == "grow", "Growth"),
                          response_var = replace(response_var, response_var == "fec",  "Fecundity") )
write.csv(best_mod_all, "results/moving_windows/best_mod_all.csv", row.names=F)



# Best models ALL COMPRISED -----------------------------------------------------------------------------

# best models (MSE)
tiff(paste0("results/moving_windows/best_mods",interval,"_mse.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(3.5,3.2,0.5,0), mgp = c(2,0.7,0) ,
    cex.lab = 1.2)

# Select best models based on MSE
mods_mse    <- mw_summ_df %>%
                  dplyr::select(species, response_var, clim_variable, model_mse) %>%
                  unique
# Set up data frame with best models
best_mod_l  <- lapply( split(mods_mse$model_mse, as.factor(mods_mse$clim_var) ),
                       function(x) table(x) )
best_mod_df <- Reduce(function(...) bind_rows(...), best_mod_l) %>% 
                  t %>% as.data.frame %>%
                  setNames( c('Air temperature','Precipitation') ) %>%
                  tibble::add_column(model = rownames(.), .before=1)

# define which colors to plot
color_df    <-  data.frame( model      = c("ctrl1", "ctrl2", "yr1", "yr2", "yr3", "yr_bet", "yr_wgt", 
                                           "expp", "gaus", 'gev', 'simpl',
                                           "expp_n", 'gev_n', 'simpl_n'),
                            model_name = c("NULL",  "24 Months", 'year_t', 'yeart-1',  'yeart-2', 'year_beta', "year_weights",
                                           "GenNorm", "Gaus", 'GEV', 'Simplex',
                                           'Gen_Norm_nested', "GEV_nested", "Simlex_nested"), 
                            col        = c("black", "grey60", 'grey80', "grey30",  "grey10", 'yellow3', 'violetred3',
                                           "green", "red", "brown", "orange",
                                           "green4", 'brown4', 'orange4') )
color_plot  <- left_join( color_df, best_mod_df )
                  
barplot(as.matrix(color_plot[,c('Air temperature','Precipitation')]), 
        beside=T, col = color_plot$col, cex.names = 1.3,
        ylab = "Number models selected as 'best'")
abline(v=15.5, lty=2)
legend("topright", color_plot$model_name, fill = color_plot$col,
       bty = "n")

dev.off()


# best models
tiff(paste0("results/moving_windows/best_mods_elpd",interval,".tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(3.5,3.2,0.5,0.5), mgp = c(2,0.7,0), cex.lab = 1.2)

# Select best models based on MSE
mods_elpd    <- mw_summ_df %>%
                  dplyr::select(species, response_var, clim_variable, model_elpd) %>%
                  unique
# set up barplot data frame
best_mod_l  <- lapply( split(mods_elpd$model_elpd, as.factor(mods_elpd$clim_var) ),
                       function(x) table(x) )
best_mod_df <- Reduce(function(...) bind_rows(...), best_mod_l) %>%
                  t %>% as.data.frame %>%
                  tibble::add_column(model = rownames(.), .before = 1 ) %>%
                  setNames( c("model","Air temperature","Precipitation") )

# define which colors to plot (color_df defined above)
color_plot  <-  full_join( color_df, best_mod_df )

barplot(as.matrix(color_plot[,c('Air temperature','Precipitation')]), 
        beside=T, col = color_plot$col, cex.names = 1.3,
        ylab = "Number models selected as 'best'")
legend("topright", color_plot$model_name, fill = color_plot$col,
       bty = "n", cex= 1)

dev.off()


# Best models BY RESPONSE ------------------------------------------------------------------------
response      <- "fec"
mw_summ_vr_df <- subset(mw_summ_df, response_var == response) %>%
                    dplyr::select(species, response_var, clim_variable, model_elpd, model_mse) %>%
                    unique

# best models
tiff(paste0("results/moving_windows/",response,"/best_mods",interval,".tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(3.5,3.2,0.5,0.5), mgp = c(2,0.7,0), cex.lab = 1.2)

# set up barplot data frame
best_mod_l  <- lapply( split(mw_summ_vr_df$model_mse, as.factor(mw_summ_vr_df$clim_var) ),
                       function(x) table(x) )
best_mod_df <- Reduce(function(...) bind_rows(...), best_mod_l) %>% 
                  t %>% as.data.frame %>%
                  tibble::add_column(model = rownames(.), .before = 1 ) %>%
                  setNames( c("model","Air temperature","Precipitation") )

# define which colors to plot
color_plot  <-  full_join( color_df, best_mod_df )

barplot(as.matrix(color_plot[,c('Air temperature','Precipitation')]), 
        beside=T, col = color_plot$col, cex.names = 1.3,
        ylab = "Number models selected as 'best'")
legend("topright", color_plot$model_name, fill = color_plot$col,
       bty = "n")

dev.off()


# best models
tiff(paste0("results/moving_windows/",response,"/best_mods_ELPD",interval,".tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(3.5,3.2,0.5,0.5), mgp = c(2,0.7,0), cex.lab = 1.2)

# set up barplot data frame
best_mod_l  <- lapply( split(mw_summ_vr_df$model_elpd, as.factor(mw_summ_vr_df$clim_var) ),
                       function(x) table(x) )
best_mod_df <- Reduce(function(...) bind_rows(...), best_mod_l) %>%
                  t %>% as.data.frame %>%
                  tibble::add_column(model = rownames(.), .before = 1 ) %>%
                  setNames( c("model","Air temperature","Precipitation") )

# define which colors to plot
color_plot   <-  full_join( color_df, best_mod_df )

# plot_mat    <- as.matrix(best_mod_df[,c("airt")])
plot_mat    <- as.matrix(best_mod_df[,c("Air temperature","Precipitation")])
barplot(as.matrix(color_plot[,c('Air temperature','Precipitation')]), 
        beside=T, col = color_plot$col, cex.names = 1.3,
        ylab = "Number models selected as 'best'")
legend("topright", color_plot$model_name, fill = color_plot$col,
       bty = "n")

dev.off()
















# best models comparison: MSE vs. ELPD
tiff(paste0("results/moving_windows/",response,"/best_mods_MSEvsELPD",interval,".tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,1), mar = c(2.5,3.2,1,0.5), mgp = c(2,0.7,0) ,
    cex.lab = 1.2)

# data frames for plotting
mse_elpd_t <- bind_rows( table(subset(mw_summ_df, clim_var == "airt")$model_mse),
                         table(subset(mw_summ_df, clim_var == "airt")$model_elpd) ) %>% 
                t %>% as.data.frame %>%
                setNames(c("MSE","ELPD")) %>%
                as.matrix
mse_elpd_p <- bind_rows( table(subset(mw_summ_df, clim_var == "precip")$model_mse),
                         table(subset(mw_summ_df, clim_var == "precip")$model_elpd) ) %>% 
                t %>% as.data.frame %>%
                setNames(c("MSE","ELPD")) %>%
                as.matrix

barplot(mse_elpd_t, beside=T, col = c("black", "grey40", "grey", "white"), 
        main = paste0("Ait temperature (",response,")"))
legend("topright", c("NULL", "24 Mon", "Expp", "Gaus"),
       fill = c("black", "grey40", "grey", "white"),
       bty = "n")
barplot(mse_elpd_p, beside=T, col = c("black", "grey40", "grey","white"), 
        main = paste0("Precipitation (",response,")") ) 

dev.off()


# MSE by absolute and year replication
tiff(paste0("results/moving_windows/",response,"/prediction_vs_rep.tiff"),
     unit="in", width=6.3, height=3.15, res=600,compression="lzw")

par(mfrow=c(1,2), mar = c(3.5,3.2,0.5,0.5), mgp = c(2,0.7,0),
    cex.lab = 1.2)
plot(jitter(mse,3) ~ jitter(rep_n,3), pch = 16, data = mw_summ_df, col = mw_summ_df$color,
     xlab = "Number of reps (year-by-site comb.)", ylab = "Mean squared error")
plot(jitter(mse,3) ~ jitter(rep_yr,3), pch = 16, data = mw_summ_df, col = mw_summ_df$color,
     xlab = "Number of years", ylab = "Mean squared error")
# legend("topleft", c("Air temp.", "GDD", "PET", "precip"), pch = 16,
#        col = unique(mw_summ_df$color), bty = "n")
legend("topright", c("Air temp.", "precip"), pch = 16,
       col = unique(mw_summ_df$color), bty = "n")

dev.off()


# MSE by climate sampled
tiff(paste0("results/moving_windows/",response,"/best_mod_MSE_by_climate_sampled.tiff"),
     unit="in", width=3.6, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,1), mar = c(3.5,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1,
    cex.axis=1)

plot(mse ~ prop_rang,xlab = "Proportion of climate values observed",
     ylab = "Best model MSE", pch = 16, data = mw_summ_df, 
     col = mw_summ_df$color)
# legend("topleft", c("Air temperature", "Growing degree days",
#                     "Potential evapotranspiration", "Precipitation"), pch = 16,
#        col = unique(mw_summ_df$color), bty = "n")
legend("topright", c("Air temperature", "Precipitation"), pch = 16,
       col = unique(mw_summ_df$color), bty = "n")
plot(mse ~ prop_yrs, xlab = "Proportion of extreme years observed",
     ylab = "Best model MSE", pch = 16, data = mw_summ_df, 
     col = mw_summ_df$color)

dev.off()


# barplots 

# color code for all (hopefully!) barplots
col_barplot <- c("red","red","green","green","orange","orange","blue","blue")

# Best models by replication
tiff(paste0("results/moving_windows/",response,"/best_mod_replication.tiff"),
     unit="in", width=6.3, height=3.15, res=600,compression="lzw")

par(mfrow=c(1,2), mar = c(2,2.5,0.1,0.5), mgp = c(1.5,0.6,0), oma =c(0,0,0,5),
    cex.axis=0.5, cex.lab = 1)
# boxplot(rep_n ~ model_climate + clim_var, data=mw_summ_df, 
#         ylab = "N. of reps (site-by-year)", names.cex=0.5,
#         names= rep(c(expression("H"[0]),"Clim."),4), col = col_barplot )
# boxplot(rep_yr ~ model_climate + clim_var, data=mw_summ_df, 
#         ylab = "N. of years sampled",
#         names= rep(c(expression("H"[0]),"Clim."),4), col = col_barplot )
# legend(8.5,36, c("Air temp.", "GDD", "PET", "precip"), xpd=NA,
#        fill = unique(col_barplot), bty = "n")
boxplot(rep_n ~ model_climate + clim_var, data=mw_summ_df, 
        ylab = "N. of reps (site-by-year)", names.cex=0.5,
        names= rep(c(expression("H"[0]),"Clim."),2), col = col_barplot )
boxplot(rep_yr ~ model_climate + clim_var, data=mw_summ_df, 
        ylab = "N. of years sampled",
        names= rep(c(expression("H"[0]),"Clim."),2), col = col_barplot )
legend(4.5,36, c("Air temp.", "Precip."), xpd=NA,
       fill = unique(col_barplot), bty = "n")

dev.off()


# best models by climate sampled
tiff(paste0("results/moving_windows/",response,"/best_mod_by_climate_sampled.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2), mar = c(2,3.5,0.5,0.2), mgp = c(2,0.5,0), oma =c(0,0,0,0),
    cex.lab = 1.2, cex.axis = 0.7)
# boxplot(prop_rang ~ model_climate + clim_var, data = mw_summ_df,
#         names = rep(c(expression("H"[0]),"Clim."),4), ylab = "Proportion of climate observed",
#         ylab = "Proportion", col = col_barplot)
# boxplot(prop_yrs ~ model_climate + clim_var, data = mw_summ_df,
#         names = rep(c(expression("H"[0]),"Clim."),4), ylab = "Proportion of extreme years",
#         ylab = "Proportion", col = col_barplot)
# boxplot(prop_var ~ model_climate + clim_var, data = mw_summ_df,
#         names = rep(c(expression("H"[0]),"Clim."),4), ylab = "Prop. of climate var. obs. (mean)",
#         ylab = "Proportion", col = col_barplot)
# boxplot(prop_var_r ~ model_climate + clim_var, data = mw_summ_df,
#         names = rep(c(expression("H"[0]),"Clim."),4), ylab = "Prop. of climate var. obs. (median)",
#         ylab = "Proportion", col = col_barplot)
# legend("bottomright", c("Air temp.", "GDD", "PET", "precip"),
#        fill = unique(col_barplot), bty = "n")
boxplot(prop_rang ~ model_climate + clim_var, data = mw_summ_df,
        names = rep(c(expression("H"[0]),"Clim."),2), ylab = "Proportion of climate observed",
        ylab = "Proportion", col = col_barplot)
boxplot(prop_yrs ~ model_climate + clim_var, data = mw_summ_df,
        names = rep(c(expression("H"[0]),"Clim."),2), ylab = "Proportion of extreme years",
        ylab = "Proportion", col = col_barplot)
boxplot(prop_var ~ model_climate + clim_var, data = mw_summ_df,
        names = rep(c(expression("H"[0]),"Clim."),2), ylab = "Prop. of climate var. obs. (mean)",
        ylab = "Proportion", col = col_barplot)
boxplot(prop_var_r ~ model_climate + clim_var, data = mw_summ_df,
        names = rep(c(expression("H"[0]),"Clim."),2), ylab = "Prop. of climate var. obs. (median)",
        ylab = "Proportion", col = col_barplot)
legend("bottomright", c("Air temp.", "Precip."),
       fill = unique(col_barplot), bty = "n")

dev.off()


# best model by "mean climate sampled" (mean climate with respect to 50 year mean)
tiff(paste0("results/moving_windows/",response,"/best_mod_by_mean_climate_sampled.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(2,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1.2)

# boxplot(mean_dev ~ model_climate + clim_var, data = mw_summ_df,
#         names = rep(c("NULL","Climate"),4), col = col_barplot,
#         ylab = "St. Dev. [mean sampled climate - mean hist. climate]" )
# abline(h = 0, lty = 2)
# legend("topleft", c("Air temperature", "Growing degree days", 
#                     "Potential evapotranspiration", "precipitation"),
#        fill = unique(col_barplot), bty = "n")

boxplot(mean_dev ~ model_climate + clim_var, data = mw_summ_df,
        names = rep(c("NULL","Climate"),2), col = col_barplot,
        ylab = "St. Dev. [mean sampled climate - mean hist. climate]" )
abline(h = 0, lty = 2)
legend("topleft", c("Air temperature","Precipitation"),
       fill = unique(col_barplot), bty = "n")

dev.off()


# species_mse by clim_var ------------------------------------------------------------
# 
# tiff(paste0("results/moving_windows/",crossval_type,"/plots/species_mse by clim_var.tiff"),
#      unit="in", width=6.3, height=6.3, res=600,compression="lzw")
# 
# # isolate only temperature and precipitation
# t_prec    <- subset(mw_summ_df, clim_var != "pet") %>%
#                 dplyr::select(clim_var, species, mse) %>%
#                 spread(clim_var,mse) %>%
#                 t
# mse       <- t_prec[-1,] %>% as.data.frame %>% setNames( t_prec[1,] )
# mse[]     <- lapply(mse[], function(x) x %<>% as.numeric(x) )
# # scale mses
# mse       <- sweep(mse, 2, STATS = colSums(mse), FUN = "/")
# mse       <- tibble::add_column(mse, x = c(1,2), .before = 1 )
# 
# # Color code - Dalgleish_spp vs. other species --------------------------------------------
# categoty  <- dplyr::select(mw_summ_df, species, Ecoregion, DicotMonoc, Class) %>%
#                 unique %>%
#                 right_join( data.frame(species = names(mse[,-1])) )
# colors    <- rep("black", (ncol(mse)-1) )
# col_dalgl <- replace( colors, names(mse[,-1]) %in% Dalgleish_spp, "red" )
# 
# par( mfrow = c(2,2), mar=c(2,3,1.5,0.1), mgp = c(2,0.7,0) )
# matplot(mse$x, mse[,-1], type = "l", main = "Dalgleish spp.",
#         ylab = "Scaled mean squared error", xlab = "",
#         col = col_dalgl, lty = 1, lwd = 1.5,
#         xaxt="n", xlim = c(0.9, 2.1))
# axis(1, at=c(1,2), labels=c("Mean temp.", "Precip."))
# legend("topright", legend = c("Other","Dalgleish"),
#        lty = 1, col = c("black","red"), bty = "n")
# 
# matplot(mse$x, mse[,-1], type = "l",  main = "Ecoregion",
#         ylab = "Scaled mean squared error", xlab = "",
#         col = categoty$Ecoregion, lty = 1, lwd = 1.5,
#         xaxt="n", xlim = c(0.9, 2.1))
# axis(1, at=c(1,2), labels=c("Mean temp.", "Precip."))
# legend("topright", legend = unique(categoty$Ecoregion),
#        lty = 1, col = 1:length(unique(categoty$Ecoregion)), bty = "n")
# 
# matplot(mse$x, mse[,-1], type = "l", main = "Dicot/Monocot",
#         ylab = "Scaled mean squared error", xlab = "",
#         col = categoty$DicotMonoc, lty = 1, lwd = 1.5,
#         xaxt="n", xlim = c(0.9, 2.1))
# axis(1, at=c(1,2), labels=c("Mean temp.", "Precip."))
# legend("topright", legend = unique(categoty$DicotMonoc),
#        lty = 1, col = 1:length(unique(categoty$DicotMonoc)), bty = "n")
# 
# matplot(mse$x, mse[,-1], type = "l", main = "Class",
#         ylab = "Scaled mean squared error", xlab = "",
#         col = categoty$Class, lty = 1, lwd = 1.5,
#         xaxt="n", xlim = c(0.9, 2.1))
# axis(1, at=c(1,2), labels=c("Mean temp.", "Precip."))
# legend("topright", legend = unique(categoty$Class),
#        lty = 1, col =1:length(unique(categoty$Class)), bty = "n")
# 
# dev.off()


# summary SPLINE plots ----------------------------------------------------------------------------------

# best models
tiff(paste0("results/splines/plots/best_mods.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(2.5,2.5,0.5,0.5), mgp = c(1.5,0.7,0) ,
    cex.lab = 1.2)

best_mod_l  <- lapply( split(spline_summ_df$model_mse, as.factor(spline_summ_df$clim_var) ),
                       function(x) table(x) )
best_mod_df <- Reduce(function(...) bind_rows(...), best_mod_l) %>% t 
colnames(best_mod_df) <- c("Air temp.","GDD","PET","Precipitation")
best_mod_df <- best_mod_df[c(2,1,3),]
barplot(best_mod_df, beside=T, col = c("black", "grey50", "grey"))
legend("topleft", c("NULL", "24 Mon", "spline"),
       fill = c("black", "grey50", "grey"),
       bty = "n")

dev.off()


# MSE by absolute and year replication
tiff(paste0("results/splines/plots/prediction_vs_rep.tiff"),
     unit="in", width=6.3, height=3.15, res=600,compression="lzw")

par(mfrow=c(1,2), mar = c(3.5,3.2,0.5,0.5), mgp = c(2,0.7,0) ,
    cex.lab = 1.2)
plot(jitter(mse,3) ~ jitter(rep_n,3), pch = 16, data = spline_summ_df, col = spline_summ_df$color,
     xlab = "Number of reps (year-by-site comb.)", ylab = "Mean squared error")
plot(jitter(mse,3) ~ jitter(rep_yr,3), pch = 16, data = spline_summ_df, col = spline_summ_df$color,
     xlab = "Number of years", ylab = "Mean squared error")
legend("topleft", c("Air temp.", "GDD", "PET", "precip"), pch = 16,
       col = unique(spline_summ_df$color), bty = "n")

dev.off()


# MSE by climate sampled
tiff(paste0("results/splines/plots/best_mod_MSE_by_climate_sampled.tiff"),
     unit="in", width=3.6, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,1), mar = c(3.5,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1,
    cex.axis=1)

plot(mse ~ prop_rang,xlab = "Proportion of climate values observed",
     ylab = "Best model MSE", pch = 16, data = spline_summ_df, 
     col = spline_summ_df$color)
legend("topleft", c("Air temperature", "Growing degree days",
                    "Potential evapotranspiration", "Precipitation"), pch = 16,
       col = unique(spline_summ_df$color), bty = "n")
plot(mse ~ prop_yrs, xlab = "Proportion of extreme years observed",
     ylab = "Best model MSE", pch = 16, data = spline_summ_df, 
     col = spline_summ_df$color)

dev.off()


# barplots 

# Best models by replication
tiff(paste0("results/splines/plots/best_mod_replication.tiff"),
     unit="in", width=6.3, height=3.15, res=600,compression="lzw")

par(mfrow=c(1,2), mar = c(2,2,0.1,0.5), mgp = c(1.1,0.5,0), oma =c(0,0,0,4.5),
    cex.axis=0.5, cex.lab = 1)
boxplot(rep_n ~ model_climate + clim_var, data=spline_summ_df, 
        ylab = "N. of reps (site-by-year)", names.cex=0.5,
        names= rep(c(expression("H"[0]),"Clim."),4), col = col_barplot )
boxplot(rep_yr ~ model_climate + clim_var, data=spline_summ_df, 
        ylab = "N. of years sampled",
        names= rep(c(expression("H"[0]),"Clim."),4), col = col_barplot )
legend(8.5,36, c("Air temp.", "GDD", "PET", "precip"), xpd=NA,
       fill = unique(col_barplot), bty = "n")

dev.off()


# best models by climate sampled
tiff(paste0("results/splines/plots/best_mod_by_climate_sampled.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2), mar = c(2,3.5,0.5,0.2), mgp = c(2,0.5,0), oma = c(0,0,0,0),
    cex.lab = 1.2, cex.axis = 0.7)

boxplot(prop_rang ~ model_climate + clim_var, data = spline_summ_df,
        names = rep(c(expression("H"[0]), "Clim."), 4), ylab = "Proportion of climate observed",
        ylab = "Proportion", col = col_barplot)
boxplot(prop_yrs ~ model_climate + clim_var, data = spline_summ_df,
        names = rep(c(expression("H"[0]), "Clim."), 4), ylab = "Proportion of extreme years",
        ylab = "Proportion", col = col_barplot)
legend("bottomright", c("Air temp.", "GDD","PET", "precip"),
       fill = unique(col_barplot), bty = "n")
boxplot(prop_var ~ model_climate + clim_var, data = spline_summ_df,
        names = rep(c(expression("H"[0]), "Clim."), 4), ylab = "Prop. of climate var. obs. (mean)",
        ylab = "Proportion", col = col_barplot)
boxplot(prop_var_r ~ model_climate + clim_var, data = spline_summ_df,
        names = rep(c(expression("H"[0]), "Clim."), 4), ylab = "Prop. of climate var. obs. (median)",
        ylab = "Proportion", col = col_barplot)

dev.off()


# best model by "mean climate sampled" (mean climate with respect to 50 year mean)
tiff(paste0("results/splines/plots/best_mod_by_mean_climate_sampled.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(2,3.5,0.5,0.2), mgp = c(2,0.7,0), oma = c(0,0,0,0), cex.lab = 1.2)

boxplot(mean_dev ~ model_climate + clim_var, data = spline_summ_df,
        names = rep(c("NULL","Climate"),4), col = col_barplot,
        ylab = "St. Dev. [mean sampled climate - mean hist. climate]" )
abline(h = 0, lty = 2)
legend("topleft", c("Air temperature", "Growing degree days","Potential evapotranspiration", "precipitation"),
       fill = unique(col_barplot), bty = "n")

dev.off()



# Models on climate/no -------------------------------------------------------------
model_climate_mods <- list(
  model_climate ~ rep_n + clim_var,
  model_climate ~ rep_yr + clim_var,
  model_climate ~ prop_rang + clim_var,
  model_climate ~ prop_yrs + clim_var,
  model_climate ~ prop_var + clim_var,
  model_climate ~ prop_var_r + clim_var,
  model_climate ~ mean_dev + clim_var,
  model_climate ~ mean_clim + clim_var,
  model_climate ~ Ecoregion + clim_var,
  model_climate ~ DicotMonoc + clim_var,
  model_climate ~ Class + clim_var,
  model_climate ~ OrganismType + clim_var
)

# fit models
mw_mods   <- lapply(model_climate_mods, function(x) glm(x, family = "binomial", data = mw_summ_df)) %>%
                    setNames(c("sample_size", "sampled_years", "prop_rang", "prop_yrs",
                               "mod_prop_var", "mod_prop_var_r", "mod_prop_mean","mod_mean_clim",
                               "ecoregion", "dicot_mono","class","organism_type"))

# summarise model results
res_summary_mod <- lapply(mw_mods, summary)


subset(mw_summ_df, model_mse != "ctrl1")$species %>% unique

# only interesting "significant" relationship
glm(model_climate ~ mean_dev, family = "binomial", 
    data = subset(mw_summ_df, clim_var == "pet") ) %>% summary


# spline models -----------------------------------------------------------------------
# fit models
sp_mods <- lapply(model_climate_mods[1:8], function(x) glm(x, family = "binomial", data = spline_summ_df)) %>%
              setNames(c("sample_size", "sampled_years", "prop_rang", "prop_yrs",
                         "mod_prop_var", "mod_prop_var_r", "mod_prop_mean","mod_mean_clim"))


lapply(sp_mods, summary)

# compare with splines --------------------------------------------------
compare_df <- merge(obs_clim_rng, mod_splin) %>%
                dplyr::select(species, model_climate, model_climate_spline)

mod <- glm(model_climate_spline ~ model_climate, family= "binomial", data = compare_df)
summary(mod)
