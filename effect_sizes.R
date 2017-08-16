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


# summarize moving windows results by climate variable
summ_by_clim <- function(clim_var){
  
  # read lambda/clim data ---------------------------------------------------------------------------------
  lam     <- read.csv("lambdas_6tr.csv", stringsAsFactors = F)
  m_info  <- read.csv("MatrixEndMonth_information.csv", stringsAsFactors = F)
  clim_fc <- data.table::fread(paste0(clim_var,"_fc_demo.csv"),  stringsAsFactors = F)
  clim_35 <- read.csv( paste0("monthly_",clim_var,"_Dalgleish.csv") )
  clim    <- list(clim_fc, clim_35)
  
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
  clim_obs_wrapper <- function(spp_name){
    
    # lambda data
    spp_lambdas   <- format_species(spp_name, lambdas)
    
    # climate data
    clim_separate <- clim_list(spp_name, clim, clim_var, spp_lambdas)
    
    # test
    expect_equal(length(spp_lambdas), length(clim_separate) )
    
    # climate ranges
    clim_rng      <- Map(observed_clim_range, clim_separate, spp_lambdas, spp_name)
    return(clim_rng)
    
  }
  # information on observed climatic ranges
  clim_rng_l_l <- lapply(spp_list, clim_obs_wrapper)
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


# plots -----------------------------------------------------------------------------------------

# climate variable and best model
tiff(paste0("results/moving_windows/loyo/plots/es/es_by_model_climVar.tiff"),
     unit="in", width=3.15, height=6.3, res=400, compression="lzw")

par(mfrow=c(2,1), mar = c(2.5,3,0.1,0.1) , mgp=c(1.8,0.7,0))
boxplot(beta ~ clim_var, ylab = "beta", data = eff_size_df) ; abline(h=0,lty=2)
boxplot(beta ~ model, ylab = "beta", data = eff_size_df) ; abline(h=0,lty=2)

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
