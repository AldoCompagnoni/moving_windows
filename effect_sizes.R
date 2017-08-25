# first try to automate graphs 
setwd("C:/cloud/MEGA/projects/sAPROPOS/")
source("C:/CODE/moving_windows/format_data.R")
library(dplyr)
library(tidyr)
library(testthat)
options(stringsAsFactors = F )

m_back        <- 24
clim_var      <- c("airt","pet","precip")
wds           <- paste0("results/moving_windows/loyo/summaries/",clim_var,"/") %>% as.list

# file list
raw_file_vec  <- lapply(wds, list.files)

# create each files' path 
path_create <- function(x,y){
  post_files <- Filter(function(x) grepl("posterior",x), x)
  paste0(y,post_files)
}

# paths of files containing posteriors
post_paths_l  <- Map(path_create, raw_file_vec, wds)
post_paths    <- Reduce(function(...) c(...), post_paths_l) %>% as.list


# calculate posterior means -------------------------------------------------------------------------
post_means <- function(x){
  
  # extract species name
  spp_n     <- gsub("results/moving_windows/loyo/summaries/[a-z]{3,6}/posterior_","",x)
  spp_name  <- gsub(".csv","",spp_n)
  
  # extract climate variable name
  # substitute 
  # 1 everything after "/posterior_"
  # 2 everything before "summaries/"
  clim_var  <- gsub("^(.*?)summaries/|/posterior_?.*","",x)
  
  posterior <- data.table::fread( paste0(x) ) %>%
                  group_by( model ) %>%
                  summarise_all( mean ) %>%
                  ungroup %>%
                  dplyr::select( -grep("log_lik_",names(.)) ) %>%
                  as.data.frame %>%
                  mutate(species = spp_name,
                         clim_var = clim_var )
  
}
means_l   <- lapply(post_paths, post_means)
means_df  <- Reduce(function(...) bind_rows(...), means_l)


# calculate posterior CI
post_ci <- function(x){
  
  # extract species name
  spp_n     <- gsub("results/moving_windows/loyo/summaries/[a-z]{3,6}/posterior_","",x)
  spp_name  <- gsub(".csv","",spp_n)
  
  # extract climate variable name
  # substitute 
  # 1 everything after "/posterior_"
  # 2 everything before "summaries/"
  clim_var  <- gsub("^(.*?)summaries/|/posterior_?.*","",x)
  
  posterior <- data.table::fread( paste0(x) ) %>%
                  dplyr::select( -grep("log_lik_",names(.)) ) %>%
                  subset( model %in% c("gaus","expp") ) %>%
                  mutate( species = spp_name,
                          clim_var = clim_var )
  
  return(posterior)
  
}
ci_l      <- lapply(post_paths, post_ci)
ci_df     <- Reduce(function(...) bind_rows(...), ci_l)



# summarize moving windows results by climate variable
summ_by_clim <- function(clim_var){
  
  # read lambda/clim data ---------------------------------------------------------------------------------
  lam     <- read.csv("lambdas_6tr.csv", stringsAsFactors = F)
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
                  subset( !is.na(MatrixEndMonth) )
  
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
                  dplyr::select( c("SpeciesAuthor","Ecoregion", "DicotMonoc", "Class") ) %>%
                  rename( species = SpeciesAuthor) %>%
                  unique
  
  # summary info ----------------------------------------------------------------------------
  res_folder<- paste0("results/moving_windows/loyo/summaries/",clim_var) 
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
  mod_splin <- read.csv( paste0("results/splines/summaries/spline_",clim_var,"24_summaries.csv") ) %>%
                  dplyr::select(species,dev0,dev1) %>%
                  unique %>%
                  mutate( model_climate = 0 ) %>%
                  mutate( model_climate = replace( model_climate, which(dev0 > dev1), 1) ) %>%
                  setNames( c("species", "dev0_spline", "dev1spline", "model_climate_spline") )
  
  # best model by measure of fit (MOF)
  best_mod_by_mof <- function(mof){
    
    lapply(mod_summ, function(x) arrange_(x, mof)) %>%
      lapply(function(x) x[1,"model"]) %>%
      unlist %>%
      table
    
  }
  best_mods <- lapply(c("mse","waic"), best_mod_by_mof) %>%
                  setNames(c("mse","waic"))
                
  # predictive accuracy by measure of fit
  pred_acc_by_mof <- function(mof){
  
    tmp <- lapply(mod_summ, function(x) arrange_(x, mof)) %>%
              lapply(function(x) x[1,c(mof,"model")]) 
    tmp <- Map(function(x,y) tibble::add_column(x, species = y, .before = 1), 
               tmp, names(tmp) )
    
     Reduce(function(...) rbind(...), tmp) %>%
        inner_join( replic ) %>%
        setNames( c("species", mof , paste0("model_", mof), paste0("rep_n_", mof)) )
  
  }
  pred_acc_l  <- Map(pred_acc_by_mof, c("mse", "looic") )
  
  # replication in years
  rep_yr <- function(spp_name){
    
    spp_df    <- subset(lambdas, SpeciesAuthor == spp_name)
    spp_yr    <- spp_df$MatrixEndYear %>% unique %>% length
    
    data.frame( species = spp_name,
                rep_yr  = spp_yr)
    
  }
  spp_rep_yr_l  <- lapply(unique(lambdas$SpeciesAuthor), rep_yr) 
  spp_rep_yr    <- Reduce(function(...) rbind(...), spp_rep_yr_l)
  
  # set up categories for summary graphs
  mod_climate <- Reduce(function(...) merge(...), pred_acc_l) %>%
                    rename( rep_n = rep_n_mse) %>%
                    left_join( spp_rep_yr ) %>%
                    inner_join( categ ) %>%
                    mutate( model_climate = sub("ctrl1", "0", model_mse)) %>%
                    mutate( model_climate = sub("ctrl2|expp|gaus", "1", model_climate)) %>%
                    mutate( model_climate = as.numeric(model_climate)) %>%
                    mutate( Ecoregion = as.factor(Ecoregion)) %>%
                    mutate( DicotMonoc = as.factor(DicotMonoc)) %>%
                    mutate( Class = as.factor(Class))
  
  # best model by climate sampled --------------------------------------------------------------
  
  # climate info 
  spp_list      <- lambdas$SpeciesAuthor %>% 
                      unique %>% 
                      sort %>% 
                      .[-12] # what happened with Daphne_rodriguezii?!?
  
  # create proportion of obeserved climates
  clim_obs_wrapper <- function(spp_name, clim_var){
    
    # lambda data
    spp_lambdas   <- format_species(spp_name, lambdas)
    
    # climate data
    clim_separate <- clim_list(spp_name, clim, spp_lambdas)
    
    # test
    expect_equal(length(spp_lambdas), length(clim_separate) )
    
    # climate ranges
    clim_rng      <- Map(observed_clim_range, clim_separate, spp_lambdas, spp_name, clim_var)
    return(clim_rng)
    
  }
  # information on observed climatic ranges
  clim_rng_l_l <- lapply(spp_list, clim_obs_wrapper, clim_var)
  clim_rng_l   <- lapply(clim_rng_l_l, function(x) Reduce(function(...) rbind(...), x))
  clim_rng_df  <- Reduce(function(...) rbind(...), clim_rng_l)
  # means of the propotion of observed climatic variability *across sites*
  clim_rng_m   <- clim_rng_df %>% group_by(species) %>% summarise_all( mean )
  obs_clim_rng <- merge(mod_climate, clim_rng_m) %>%
                    mutate( clim_var = clim_var )
  
  return( obs_clim_rng )
  
}
# put it all together
mw_summ_l   <- lapply(c("precip", "pet", "airt"), summ_by_clim )
# Also: color code the data frame 
mw_summ_df  <- Reduce(function(...) rbind(...), mw_summ_l) %>%
                  merge( data.frame( clim_var = c("precip", "pet", "airt"),
                                     color    = c("blue", "orange", "red") ) ) %>%
                  mutate(method = "mov_win")
                  

# effect sizes ----------------------------------------------------------------------------------
eff_size_df <- mw_summ_df %>% 
                  dplyr::select( -one_of( c("looic","model_looic","model_looic","rep_n_looic") ) ) %>%
                  rename( model = model_mse ) %>%
                  left_join(means_df) %>%
                  subset(model != "ctrl1" )

# effect sizes - control *beta* of ctrl2 when null model (ctrl1) is the best model
means_0     <- means_df %>% 
                  subset( model == "ctrl2" ) %>%
                  mutate( model =  replace(model, model == "ctrl2", "ctrl1") ) 
eff_size_0  <- mw_summ_df %>% 
                  dplyr::select( -one_of( c("looic","model_looic","model_looic","rep_n_looic") ) ) %>%
                  rename( model = model_mse ) %>%
                  subset(model == "ctrl1" ) %>%
                  left_join(means_0) 

# all effect sizes
all_e_s_df  <- bind_rows(eff_size_df, eff_size_0) %>% arrange(clim_var, species)


# plots -----------------------------------------------------------------------------------------

# climate variable and best model
tiff(paste0("results/moving_windows/loyo/plots/es/es_by_model_climVar.tiff"),
     unit="in", width=3.15, height=6.3, res=400, compression="lzw")

par(mfrow=c(2,1), mar = c(2.2,3,0.1,0.1) , mgp=c(1.8,0.7,0))
boxplot(beta ~ clim_var, ylab = "Effect of climate (beta)", data = eff_size_df) ; abline(h=0,lty=2)
boxplot(beta ~ model, ylab = "Effect of climate (beta)", data = all_e_s_df) ; abline(h=0,lty=2)

dev.off()


# Categories
tiff(paste0("results/moving_windows/loyo/plots/es/es_by_category.tiff"),
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par(mfrow=c(2,2), mar = c(2.5,3,0.1,0.1) , mgp=c(1.8,0.7,0))
boxplot(beta ~ Ecoregion, ylab = "beta", data = eff_size_df) ; abline(h=0,lty=2)
boxplot(beta ~ DicotMonoc, ylab = "beta", data = eff_size_df) ; abline(h=0,lty=2)
boxplot(beta ~ Class, ylab = "beta", data = eff_size_df) ; abline(h=0,lty=2)

dev.off()

# replication
tiff(paste0("results/moving_windows/loyo/plots/es/es_by_replication.tiff"),
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par(mfrow=c(2,2), mar = c(2.5,3,0.1,0.1) , mgp=c(1.8,0.7,0))
plot(beta ~ mse, pch = 16,ylab="beta", col = color, data = eff_size_df)
plot(beta ~ rep_n, pch = 16, ylab="beta",col = color, data = eff_size_df)
plot(beta ~ rep_yr, pch = 16, ylab="beta",col = color, data = eff_size_df)
dev.off()

# proportion of climate sampled
tiff(paste0("results/moving_windows/loyo/plots/es/es_by_climateSampled.tiff"),
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par(mfrow=c(2,2), mar = c(2.5,3,0.1,0.1), mgp=c(1.8,0.7,0))
plot(beta ~ prop_rang, pch = 16, ylab="beta",col = color, data = eff_size_df)
plot(beta ~ prop_yrs, pch = 16, ylab="beta",col = color, data = eff_size_df)
plot(beta ~ prop_var, pch = 16, ylab="beta",col = color, data = eff_size_df)
plot(beta ~ prop_var_r, pch = 16, ylab="beta",col = color, data = eff_size_df)
dev.off()

# mean climate sampled
tiff(paste0("results/moving_windows/loyo/plots/es/es_by_meanClimate.tiff"),
     unit="in", width=3.15, height=6.3, res=400, compression="lzw")

par(mfrow=c(2,1), mar = c(3,3,0.1,0.1), mgp=c(1.8,0.7,0))
plot(beta ~ mean_dev, pch = 16, ylab="beta", col = color, data = eff_size_df)
plot(beta ~ mean_clim, pch = 16,ylab="beta", col = color, data = eff_size_df)
dev.off()





# isolate only temperature and precipitation -------------------------------------------------
t_prec    <- subset(eff_size_df, clim_var != "pet") %>%
                dplyr::select(clim_var, species, beta) %>%
                spread(clim_var,beta) %>% 
                t 
beta      <- t_prec[-1,] %>% as.data.frame %>% setNames( t_prec[1,] )
beta[]    <- lapply(beta[], function(x) x %<>% as.numeric(x) )
# scale mses
beta      <- tibble::add_column(beta, x = c(1,2), .before = 1 )


# Color code - Dalgleish_spp vs. other species --------------------------------------------
categoty  <- dplyr::select(mw_summ_df, species, Ecoregion, DicotMonoc, Class) %>% 
                unique %>%
                right_join( data.frame(species = names(beta[,-1])) )
colors    <- rep("black", (ncol(beta)-1) )
col_dalgl <- replace( colors, names(beta[,-1]) %in% Dalgleish_spp, "red" )


par( mfrow = c(2,2), mar=c(2,3,1.5,0.1), mgp = c(2,0.7,0) )
matplot(beta$x, beta[,-1], type = "l", main = "Dalgleish spp.",
        ylab = "Beta", xlab = "",
        col = col_dalgl, lty = 1, lwd = 1.5,
        xaxt="n", xlim = c(0.9, 2.1))
axis(1, at=c(1,2), labels=c("Mean temp.", "Precip.")) 
legend("topright", legend = c("Other","Dalgleish"),
       lty = 1, col = c("black","red"), bty = "n")

matplot(beta$x, beta[,-1], type = "l",  main = "Ecoregion",
        ylab = "Beta", xlab = "",
        col = categoty$Ecoregion, lty = 1, lwd = 1.5,
        xaxt="n", xlim = c(0.9, 2.1))
axis(1, at=c(1,2), labels=c("Mean temp.", "Precip.")) 
legend("topright", legend = unique(categoty$Ecoregion),
       lty = 1, col = 1:length(unique(categoty$Ecoregion)), bty = "n")

matplot(beta$x, beta[,-1], type = "l", main = "Dicot/Monocot",
        ylab = "Beta", xlab = "",
        col = categoty$DicotMonoc, lty = 1, lwd = 1.5,
        xaxt="n", xlim = c(0.9, 2.1))
axis(1, at=c(1,2), labels=c("Mean temp.", "Precip.")) 
legend("topright", legend = unique(categoty$DicotMonoc),
       lty = 1, col = 1:length(unique(categoty$DicotMonoc)), bty = "n")

matplot(beta$x, beta[,-1], type = "l", main = "Class",
        ylab = "Beta", xlab = "",
        col = categoty$Class, lty = 1, lwd = 1.5,
        xaxt="n", xlim = c(0.9, 2.1))
axis(1, at=c(1,2), labels=c("Mean temp.", "Precip.")) 
legend("topright", legend = unique(categoty$Class),
       lty = 1, col =1:length(unique(categoty$Class)), bty = "n")




# Credible Intervals of location and shape ----------------------------------------------------------------

# upper and lower Credible Intervals
ci_05     <- ci_df %>%
                group_by( model, species, clim_var ) %>%
                summarise_all( quantile_05 ) %>%
                setNames( c("model", "species", "clim_var",
                            paste0( c("sens_mu", "sens_sd", "alpha", "beta", "y_sd", "lp__"),"95") ) )
ci_95     <- ci_df %>%
                group_by( model, species, clim_var ) %>%
                summarise_all( quantile_95 ) %>%
                setNames( c("model", "species", "clim_var",
                            paste0( c("sens_mu", "sens_sd", "alpha", "beta", "y_sd", "lp__"),"05") ) )
up_lo_ci  <- full_join(ci_05,ci_95) 



# draw bargraph and colors
bars_and_ci <- function(variable, climate_var, models){
  
  # set 'base' colors
  main_col  <- data.frame( clim_var = c("airt","pet","precip"),
                           col      = c("red","orange","blue") ) %>%
                  subset( clim_var == climate_var ) %>% .$col
  
  # set up data frame of colors
  col_df    <- mw_summ_df %>% 
                  dplyr::select(species, clim_var, model_mse) %>%
                  arrange(clim_var, species) %>%
                  subset( clim_var == climate_var ) %>%
                  mutate( col = main_col ) %>%
                  mutate( col = replace(col, model_mse == models,  paste0(main_col,"4") ) )
  
  if( variable == "sens_mu"){
    ylim_vec=c(0,24)
  # if measureing
  }else{
    if(models == "gaus"){
      ylim_vec=c(0,45)
    }else{
      ylim_vec=c(0,12)
    }
  }
    
  # plot interquartile range
  boxplot( formula( paste0(variable," ~ species") ), names = NA, ylab = "", 
          outline = F, ylim = ylim_vec, whisklty = 0, staplelty = 0, xaxt='n', 
          col = col_df$col, 
          data = subset(ci_df, model == models & clim_var == climate_var))
  
  # upper and lower credible intervals
  up_lo_df  <- up_lo_ci %>%
                  subset(model == models & clim_var == climate_var) %>%
                  .[,c("species",paste0(variable,c("05","95")))] %>%
                  as.data.frame
  # species list
  spp_ci    <- up_lo_df$species %>% unique
  
  # plot upper and lower CI - by species
  plot_ci   <- function(ii){  
  
    segOff  <- 0.5
    # upper ci
    segments(x0=ii-segOff,y0 = up_lo_df[ii,2], x1=ii+segOff, y1= up_lo_df[ii,2], lty=1, lwd=1.5, col="black")
    # lower ci
    segments(x0=ii-segOff,y0 = up_lo_df[ii,3], x1=ii+segOff, y1= up_lo_df[ii,3], lty=1, lwd=1.5, col="black")
    # green vertical line
    segments(x0=ii,y0 = up_lo_df[ii,3], x1=ii, y1=up_lo_df[ii,2], lty=3, lwd=1, col="black")
      
  }
  lapply(1:length(spp_ci), plot_ci)
    
}


# "center of sensitivity" for both 
tiff("results/moving_windows/loyo/plots/es/central_tendency_ci.tiff",unit="in",
     width=6.3,height=4.2,res=600,compression="lzw")

par(mfrow = c(2,3), mar = c(2,2,0.1,0.1), mgp = c(2,0.7,0), oma = c(0,1,0,0) )
# gaussian models
bars_and_ci("sens_mu", "airt","gaus")
bars_and_ci("sens_mu", "pet","gaus")
bars_and_ci("sens_mu", "precip","gaus")

# power exponential
bars_and_ci("sens_mu", "airt","expp")
bars_and_ci("sens_mu", "pet","expp")
bars_and_ci("sens_mu", "precip","expp")

# x-axis
mtext("Species", side = 1, line = 0.7, at = -24)

# y-axes
mtext("Gaussian: mean (month)",              line = 32.9, at = 42, cex = 0.8, side = 2)
mtext("Power exponential: location (month)", line = 32.9, at = 11, cex = 0.8, side = 2)

# plot data
mtext("Air Temperature", side = 1, line = -15.3, at = -67)
mtext("PET",             side = 1, line = -15.3, at = -26)
mtext("Precipitation",   side = 1, line = -15.3, at = 18)

dev.off()


# standard deviation; scale
tiff("results/moving_windows/loyo/plots/es/dispersion_ci.tiff",unit="in",
     width=6.3,height=4.2,res=600,compression="lzw")

par(mfrow = c(2,3), mar = c(2,2,0.1,0.1), mgp = c(2,0.7,0), oma = c(0,2,0,0) )

# gaussian models
bars_and_ci("sens_sd", "airt","gaus")
bars_and_ci("sens_sd", "pet","gaus")
bars_and_ci("sens_sd", "precip","gaus")

# power exponential
bars_and_ci("sens_sd", "airt","expp")
bars_and_ci("sens_sd", "pet","expp")
bars_and_ci("sens_sd", "precip","expp")

# x-axis
mtext("Species", side = 1, line = 0.7, at = -24)

# y-axes
mtext("Power exponential: location (month)",line = 32.5, at = 11, cex = 0.8, side = 2 )
mtext("Gaussian: location (month)",         line = 33.2, at = 42, cex = 0.8, side = 2 )

# plot data
mtext("Air Temperature", side = 1, line = -15.3, at = -67)
mtext("PET",             side = 1, line = -15.3, at = -26)
mtext("Precipitation",   side = 1, line = -15.3, at = 18)

dev.off()


# standard deviation; scale
tiff("results/moving_windows/loyo/plots/es/legend.tiff",unit="in",
     width=6.3,height=6.3,res=600,compression="lzw")

par(mfrow = c(1,1), mar = c(2,2,0.1,0.1), mgp = c(2,0.7,0), oma = c(0,2,0,0) )
plot(c(1:10), type="n", xaxt = "n", yaxt = "n", bty='n')
legend(0.2,10,
       c("Air temperature - not best model",
         "Air temperature - BEST model!",
         "PET - not best model",
         "PET - BEST model!",
         "Precipitation - not best model",
         "Precipitation - BEST model!"),
       fill = c("red","red4","orange","orange4","blue","blue4"),
       bty="n", cex = 2)

dev.off()
