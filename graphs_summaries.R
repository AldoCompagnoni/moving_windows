# first try to automate graphs 
setwd("C:/cloud/MEGA/projects/sAPROPOS/")
source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(magrittr)
library(testthat)
options(stringsAsFactors = F )

crossval_type <- "loyo"
m_back        <- 24


# summarize moving windows results by climate variable
summ_by_clim <- function(clim_var){
  
  # read lambda/clim data ---------------------------------------------------------------------------------
  lam     <- read.csv("lambdas_6tr.csv", stringsAsFactors = F)
  m_info  <- read.csv("MatrixEndMonth_information.csv", stringsAsFactors = F)
  clim_fc   <- data.table::fread(paste0(clim_var,"_fc_demo.csv"),  stringsAsFactors = F)
  clim_35   <- read.csv( paste0("monthly_",clim_var,"_Dalgleish.csv") )
  clim      <- list(clim_fc, clim_35)
  
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
  res_folder<- paste0("results/moving_windows/",crossval_type,"/summaries/",clim_var) 
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
spline_summ_l   <- lapply( c("airt","pet","precip"), read_spline_summ, 12)
spline_summ_tmp <- Reduce( function(...) rbind(...), spline_summ_l ) %>%
                      merge( data.frame( clim_var = c("precip", "pet", "airt"),
                                         color    = c("blue", "orange", "red") ) ) %>%
                      mutate( method = "spline") %>%
                      rowwise %>%
                      mutate(mse = min(dev_NULL,dev_lm,dev_spline) )
                      
# add information on replication, and climate sampled
rep_clim_info   <- dplyr::select(mw_summ_df, species, clim_var, rep_n, rep_yr,
                                 prop_rang, prop_yrs,prop_var,prop_var_r, mean_dev, mean_clim ) %>%
                      unique
spline_summ_df  <- left_join( spline_summ_tmp, rep_clim_info )


# summary MOVING WINDOWS plots ----------------------------------------------------------------------------------

# best models
tiff(paste0("results/moving_windows/",crossval_type,"/plots/best_mods.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(3.5,3.2,0.5,0.5), mgp = c(2,0.7,0) ,
    cex.lab = 1.2)

best_mod_l  <- lapply( split(mw_summ_df$model_mse, as.factor(mw_summ_df$clim_var) ),
                       function(x) table(x) )
best_mod_df <- Reduce(function(...) bind_rows(...), best_mod_l) %>% t 
colnames(best_mod_df) <- c("Air temp.","PET","Precipitation")
barplot(best_mod_df, beside=T, col = c("black", "grey40", "grey", "white"))
legend("topright", c("NULL", "24 Mon", "Expp", "Gaus"),
       fill = c("black", "grey40", "grey", "white"),
       bty = "n")

dev.off()


# MSE by absolute and year replication
tiff(paste0("results/moving_windows/",crossval_type,"/plots/prediction_vs_rep.tiff"),
     unit="in", width=6.3, height=3.15, res=600,compression="lzw")

par(mfrow=c(1,2), mar = c(3.5,3.2,0.5,0.5), mgp = c(2,0.7,0) ,
    cex.lab = 1.2)
plot(jitter(mse,3) ~ jitter(rep_n,3), pch = 16, data = mw_summ_df, col = mw_summ_df$color,
     xlab = "Number of reps (year-by-site comb.)", ylab = "Mean squared error")
plot(jitter(mse,3) ~ jitter(rep_yr,3), pch = 16, data = mw_summ_df, col = mw_summ_df$color,
     xlab = "Number of years", ylab = "Mean squared error")
legend("topleft", c("Air temp.", "PET", "precip"), pch = 16,
       col = unique(mw_summ_df$color), bty = "n")

dev.off()


# MSE by climate sampled
tiff(paste0("results/moving_windows/",crossval_type,"/plots/best_mod_MSE_by_climate_sampled.tiff"),
     unit="in", width=3.6, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,1), mar = c(3.5,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1,
    cex.axis=1)

plot(mse ~ prop_rang,xlab = "Proportion of climate values observed",
     ylab = "Best model MSE", pch = 16, data = mw_summ_df, 
     col = mw_summ_df$color)
legend("topleft", c("Air temperature", "Potential evapotranspiration", "Precipitation"), pch = 16,
       col = unique(mw_summ_df$color), bty = "n")
plot(mse ~ prop_yrs, xlab = "Proportion of extreme years observed",
     ylab = "Best model MSE", pch = 16, data = mw_summ_df, 
     col = mw_summ_df$color)

dev.off()


# barplots 

# color code for all (hopefully!) barplots
col_barplot <- c("red","red","orange","orange","blue","blue")

# Best models by replication
tiff(paste0("results/moving_windows/",crossval_type,"/plots/best_mod_replication.tiff"),
     unit="in", width=6.3, height=3.15, res=600,compression="lzw")

par(mfrow=c(1,2), mar = c(2,2.5,0.1,0.5), mgp = c(1.5,0.6,0), oma =c(0,0,0,5),
    cex.axis=0.5, cex.lab = 1)
boxplot(rep_n ~ model_climate + clim_var, data=mw_summ_df, 
        ylab = "N. of reps (site-by-year)", names.cex=0.5,
        names= rep(c("NULL","Clim."),3), col = col_barplot )
boxplot(rep_yr ~ model_climate + clim_var, data=mw_summ_df, 
        ylab = "N. of years sampled",
        names= rep(c("NULL","Clim."),3), col = col_barplot )
legend(6.5,36, c("Air temp.", "PET", "precip"), xpd=NA,
       fill = unique(col_barplot), bty = "n")

dev.off()

# best models by climate sampled
tiff(paste0("results/moving_windows/",crossval_type,"/plots/best_mod_by_climate_sampled.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2), mar = c(2,3.5,0.5,0.2), mgp = c(2,0.5,0), 
    cex.lab = 1.2, cex.axis = 0.7)

boxplot(prop_rang ~ model_climate + clim_var, data = mw_summ_df,
        names = rep(c("NULL", "Climate"), 3), ylab = "Proportion of climate observed",
        ylab = "Proportion", col = col_barplot)
boxplot(prop_yrs ~ model_climate + clim_var, data = mw_summ_df,
        names = rep(c("NULL", "Climate"), 3), ylab = "Proportion of extreme years",
        ylab = "Proportion", col = col_barplot)
legend("bottomright", c("Air temp.", "PET", "precip"),
       fill = unique(col_barplot), bty = "n")
boxplot(prop_var ~ model_climate + clim_var, data = mw_summ_df,
        names = rep(c("NULL", "Climate"), 3), ylab = "Prop. of climate var. obs. (mean)",
        ylab = "Proportion", col = col_barplot)
boxplot(prop_var_r ~ model_climate + clim_var, data = mw_summ_df,
        names = rep(c("NULL", "Climate"), 3), ylab = "Prop. of climate var. obs. (median)",
        ylab = "Proportion", col = col_barplot)

dev.off()


# best model by "mean climate sampled" (mean climate with respect to 50 year mean)
tiff(paste0("results/moving_windows/",crossval_type,"/plots/best_mod_by_mean_climate_sampled.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(2,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1.2)

boxplot(mean_dev ~ model_climate + clim_var, data = mw_summ_df,
        names = rep(c("NULL","Climate"),3), col = col_barplot,
        ylab = "St. Dev. [mean sampled climate - mean hist. climate]" )
abline(h = 0, lty = 2)
legend("topleft", c("Air temperature", "Potential evapotranspiration", "precipitation"),
       fill = unique(col_barplot), bty = "n")

dev.off()


# species_mse by clim_var ------------------------------------------------------------

tiff(paste0("results/moving_windows/",crossval_type,"/plots/species_mse by clim_var.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

# isolate only temperature and precipitation 
t_prec    <- subset(mw_summ_df, clim_var != "pet") %>%
                dplyr::select(clim_var, species, mse) %>%
                spread(clim_var,mse) %>% 
                t 
mse       <- t_prec[-1,] %>% as.data.frame %>% setNames( t_prec[1,] )
mse[]     <- lapply(mse[], function(x) x %<>% as.numeric(x) )
# scale mses
mse       <- sweep(mse, 2, STATS = colSums(mse), FUN = "/")
mse       <- tibble::add_column(mse, x = c(1,2), .before = 1 )

# Color code - Dalgleish_spp vs. other species --------------------------------------------
categoty  <- dplyr::select(mw_summ_df, species, Ecoregion, DicotMonoc, Class) %>% 
                unique %>%
                right_join( data.frame(species = names(mse[,-1])) )
colors    <- rep("black", (ncol(mse)-1) )
col_dalgl <- replace( colors, names(mse[,-1]) %in% Dalgleish_spp, "red" )
  
par( mfrow = c(2,2), mar=c(2,3,1.5,0.1), mgp = c(2,0.7,0) )
matplot(mse$x, mse[,-1], type = "l", main = "Dalgleish spp.",
        ylab = "Scaled mean squared error", xlab = "",
        col = col_dalgl, lty = 1, lwd = 1.5,
        xaxt="n", xlim = c(0.9, 2.1))
axis(1, at=c(1,2), labels=c("Mean temp.", "Precip.")) 
legend("topright", legend = c("Other","Dalgleish"),
       lty = 1, col = c("black","red"), bty = "n")

matplot(mse$x, mse[,-1], type = "l",  main = "Ecoregion",
        ylab = "Scaled mean squared error", xlab = "",
        col = categoty$Ecoregion, lty = 1, lwd = 1.5,
        xaxt="n", xlim = c(0.9, 2.1))
axis(1, at=c(1,2), labels=c("Mean temp.", "Precip.")) 
legend("topright", legend = unique(categoty$Ecoregion),
       lty = 1, col = 1:length(unique(categoty$Ecoregion)), bty = "n")

matplot(mse$x, mse[,-1], type = "l", main = "Dicot/Monocot",
        ylab = "Scaled mean squared error", xlab = "",
        col = categoty$DicotMonoc, lty = 1, lwd = 1.5,
        xaxt="n", xlim = c(0.9, 2.1))
axis(1, at=c(1,2), labels=c("Mean temp.", "Precip.")) 
legend("topright", legend = unique(categoty$DicotMonoc),
       lty = 1, col = 1:length(unique(categoty$DicotMonoc)), bty = "n")

matplot(mse$x, mse[,-1], type = "l", main = "Class",
        ylab = "Scaled mean squared error", xlab = "",
        col = categoty$Class, lty = 1, lwd = 1.5,
        xaxt="n", xlim = c(0.9, 2.1))
axis(1, at=c(1,2), labels=c("Mean temp.", "Precip.")) 
legend("topright", legend = unique(categoty$Class),
       lty = 1, col =1:length(unique(categoty$Class)), bty = "n")

dev.off()


# summary SPLINE plots ----------------------------------------------------------------------------------

# best models
tiff(paste0("results/splines/plots/best_mods.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(2.5,2.5,0.5,0.5), mgp = c(1.5,0.7,0) ,
    cex.lab = 1.2)

best_mod_l  <- lapply( split(spline_summ_df$model_mse, as.factor(spline_summ_df$clim_var) ),
                       function(x) table(x) )
best_mod_df <- Reduce(function(...) bind_rows(...), best_mod_l) %>% t 
colnames(best_mod_df) <- c("Air temp.","PET","Precipitation")
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
legend("topleft", c("Air temp.", "PET", "precip"), pch = 16,
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
legend("topleft", c("Air temperature", "Potential evapotranspiration", "Precipitation"), pch = 16,
       col = unique(spline_summ_df$color), bty = "n")
plot(mse ~ prop_yrs, xlab = "Proportion of extreme years observed",
     ylab = "Best model MSE", pch = 16, data = spline_summ_df, 
     col = spline_summ_df$color)

dev.off()


# barplots 

# Best models by replication
tiff(paste0("results/splines/plots/best_mod_replication.tiff"),
     unit="in", width=6.3, height=3.15, res=600,compression="lzw")

par(mfrow=c(1,2), mar = c(2,2.5,0.1,0.5), mgp = c(1.5,0.6,0), oma =c(0,0,0,5),
    cex.axis=0.5, cex.lab = 1)
boxplot(rep_n ~ model_climate + clim_var, data=spline_summ_df, 
        ylab = "N. of reps (site-by-year)", names.cex=0.5,
        names= rep(c("NULL","Clim."),3), col = col_barplot )
boxplot(rep_yr ~ model_climate + clim_var, data=spline_summ_df, 
        ylab = "N. of years sampled",
        names= rep(c("NULL","Clim."),3), col = col_barplot )
legend(6.5,36, c("Air temp.", "PET", "precip"), xpd=NA,
       fill = unique(col_barplot), bty = "n")

dev.off()


# best models by climate sampled
tiff(paste0("results/splines/plots/best_mod_by_climate_sampled.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2), mar = c(2,3.5,0.5,0.2), mgp = c(2,0.5,0), 
    cex.lab = 1.2, cex.axis = 0.7)

boxplot(prop_rang ~ model_climate + clim_var, data = spline_summ_df,
        names = rep(c("NULL", "Climate"), 3), ylab = "Proportion of climate observed",
        ylab = "Proportion", col = col_barplot)
boxplot(prop_yrs ~ model_climate + clim_var, data = spline_summ_df,
        names = rep(c("NULL", "Climate"), 3), ylab = "Proportion of extreme years",
        ylab = "Proportion", col = col_barplot)
legend("bottomright", c("Air temp.", "PET", "precip"),
       fill = unique(col_barplot), bty = "n")
boxplot(prop_var ~ model_climate + clim_var, data = spline_summ_df,
        names = rep(c("NULL", "Climate"), 3), ylab = "Prop. of climate var. obs. (mean)",
        ylab = "Proportion", col = col_barplot)
boxplot(prop_var_r ~ model_climate + clim_var, data = spline_summ_df,
        names = rep(c("NULL", "Climate"), 3), ylab = "Prop. of climate var. obs. (median)",
        ylab = "Proportion", col = col_barplot)

dev.off()


# best model by "mean climate sampled" (mean climate with respect to 50 year mean)
tiff(paste0("results/splines/plots/best_mod_by_mean_climate_sampled.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(2,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1.2)

boxplot(mean_dev ~ model_climate + clim_var, data = spline_summ_df,
        names = rep(c("NULL","Climate"),3), col = col_barplot,
        ylab = "St. Dev. [mean sampled climate - mean hist. climate]" )
abline(h = 0, lty = 2)
legend("topleft", c("Air temperature", "Potential evapotranspiration", "precipitation"),
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
  model_climate ~ Class + clim_var
)

# fit models
mw_mods   <- lapply(model_climate_mods, function(x) glm(x, family = "binomial", data = mw_summ_df)) %>%
                    setNames(c("sample_size", "sampled_years", "prop_rang", "prop_yrs",
                               "mod_prop_var", "mod_prop_var_r", "mod_prop_mean","mod_mean_clim",
                               "ecoregion", "dicot_mono","class"))

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
