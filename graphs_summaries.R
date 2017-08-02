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
    select( c("SpeciesAuthor", var) ) %>%
    distinct %>%
    group_by_(var) %>%
    summarise(spp_n = n())

}
# replication of categories
categ_rep <- lapply(c("OrganismType","Ecoregion","Order", "Class", "DicotMonoc", "Altitude"), 
                    by_spp_rep)

# I chose to use the following categories
categ     <- lambdas %>%
                select( c("SpeciesAuthor","Ecoregion", "DicotMonoc", "Class") ) %>%
                rename( species = SpeciesAuthor) %>%
                unique


# summary info ----------------------------------------------------------------------------
res_folder<- "supercomp/res_7.26" 
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

# best model by measure of fit (MOF)
best_mod_by_mof <- function(mof){
  
  lapply(mod_summ, function(x) arrange_(x, mof)) %>%
    lapply(function(x) x[1,"model"]) %>%
    unlist %>%
    table
  
}
best_mods <- lapply(c("mse","looic"), best_mod_by_mof) %>%
                setNames(c("mse","looic"))
              
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

# set up categories for summary graphs
mod_climate <- Reduce(function(...) merge(...), pred_acc_l) %>%
                  rename( rep_n = rep_n_mse) %>%
                  inner_join( categ ) %>%
                  mutate( model_climate = sub("ctrl1", "0", model_mse)) %>%
                  mutate( model_climate = sub("ctrl2|expp", "1", model_climate)) %>%
                  mutate( model_climate = as.numeric(model_climate)) %>%
                  mutate( Ecoregion = as.factor(Ecoregion)) %>%
                  mutate( DicotMonoc = as.factor(DicotMonoc)) %>%
                  mutate( Class = as.factor(Class))

# calculate number of categories for each "best model"
best_mod_by_categ <- function(best_mod = NULL, category){
  
  if( is.null(best_mod) ){
    tmp <- mod_climate %>%
      group_by_( category ) %>%
      summarise( rep = n() ) %>%
      as.data.frame %>% t
  } else{
    tmp <- mod_climate %>%
      subset( model_climate == best_mod ) %>%
      group_by_( category ) %>%
      summarise( rep = n() ) %>%
      as.data.frame %>% t
  }
  df <- matrix(tmp[2,], nrow=1, ncol=length(tmp[2,]), byrow=T) %>%
    as.data.frame %>%
    setNames( tmp[1,] ) %>% 
    as.matrix()
  return(df)
  
}


# summary plots ----------------------------------------------------------------------------------

# best models
tiff(paste0("results/plots/best_mods.tiff"),
     unit="in", width=6.3, height=3.15, res=600,compression="lzw")

par(mfrow=c(1,2), mar = c(3.5,3.2,0.5,0.5), mgp = c(2,0.7,0) ,
    cex.lab = 1.2)

best_mods$mse   <- setNames(best_mods$mse, c("NULL", "24mon.", "expp"))
best_mods$looic <- setNames(best_mods$looic, "Gaussian")
barplot(best_mods$mse, ylab = "Best model based on MSE", cex.names = 1.2)
barplot(best_mods$looic, ylab = "Best model based on LOOIC", cex.names = 1.2)

dev.off()


# prediction vs. replication
tiff(paste0("results/plots/prediction_vs_rep.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(3.5,3.2,0.5,0.5), mgp = c(2,0.7,0) ,
    cex.lab = 1.2)
plot(mse ~ rep_n, pch = 16, data = mod_climate, 
     xlab = "Number of replicates", ylab = "Mean squared error")

dev.off()


# replication of best models
tiff(paste0("results/plots/best_mod_replication.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2), mar = c(3.5,3.5,1,0.5), mgp = c(2,0.7,0) ,
    cex.lab = 1.2)

hist(mod_climate$rep_n, ylim = c(0,8), xlim = c(0,70), xlab = "replication", main = "All data sets")
hist(subset(mod_climate, model_climate == 0)$rep_n, ylim = c(0,8), xlim = c(0,70), xlab = "replication", main = "No climate")
hist(subset(mod_climate, model_climate == 1)$rep_n, ylim = c(0,8), xlim = c(0,70), xlab = "replication", main = "Climate")

dev.off()


# best model by ecoregion
tiff(paste0("results/plots/best_mod_by_ecoregion.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2), mar = c(3,3.5,1.5,0.5), mgp = c(2,0.7,0), cex.lab = 1.2)
barplot( best_mod_by_categ(NULL,"Ecoregion"), main = "Representation across species", 
         ylim = c(0,14), col = "grey")
barplot( best_mod_by_categ(1,"Ecoregion"), main = "Climate", ylim = c(0,14), col = "grey")
barplot( best_mod_by_categ(0,"Ecoregion"), main = "NULL", ylim = c(0,14), col = "grey")

dev.off()


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
  clim_separate <- clim_list(spp_name, clim, spp_lambdas)
  
  # test
  expect_true( length(spp_lambdas) == length(clim_separate) )
  
  # climate ranges
  clim_rng      <- Map(observed_clim_range, clim_separate, spp_lambdas, spp_name)
  clim_rng
  
}
# information on observed climatic ranges
clim_rng_l_l <- lapply(spp_list, clim_obs_wrapper)
clim_rng_l   <- lapply(clim_rng_l_l, function(x) Reduce(function(...) rbind(...), x))
clim_rng_df  <- Reduce(function(...) rbind(...), clim_rng_l)
# means of the propotion of observed climatic variability *across sites*
clim_rng_m   <- clim_rng_df %>% group_by(species) %>% summarise_all( mean )



# best model by observed climate range -----------------------------------------
obs_clim_rng <- merge(mod_climate, clim_rng_m)

tiff(paste0("results/plots/best_mod_by_climate_sampled.tiff"),
     unit="in", width=3.6, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2), mar = c(2,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1.2)

boxplot(prop_rang ~ model_climate, data = obs_clim_rng,
        names = c("No Climate", "Climate"), ylab = "Proportion of climate observed",
        ylab = "Proportion")
boxplot(prop_yrs ~ model_climate, data = obs_clim_rng,
        names = c("No Climate", "Climate"), ylab = "Proportion of extreme years",
        ylab = "Proportion")
boxplot(prop_var ~ model_climate, data = obs_clim_rng,
        names = c("No Climate", "Climate"), ylab = "Prop. of climate var. obs. (mean)",
        ylab = "Proportion")
boxplot(prop_var_r ~ model_climate, data = obs_clim_rng,
        names = c("No Climate", "Climate"), ylab = "Prop. of climate var. obs. (median)",
        ylab = "Proportion")

dev.off()


tiff(paste0("results/plots/best_mod_by_mean_climate_sampled.tiff"),
     unit="in", width=3.6, height=6.3, res=600,compression="lzw")

par(mfrow=c(1,1), mar = c(2,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1.2)

boxplot(mean_dev ~ model_climate, data = obs_clim_rng,
        names = c("No Climate", "Climate"), 
        ylab = "St. Dev. of mean sampled climate from mean hist. climate" )
abline(h = 0, lty = 2)

dev.off()


tiff(paste0("results/plots/best_mod_MSE_by_climate_sampled.tiff"),
     unit="in", width=3.6, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,1), mar = c(3.5,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1.2)

plot(mse ~ prop_rang,xlab = "Proportion of climate values observed",
     ylab = "Best model MSE", pch = 16, data = obs_clim_rng, col = "black")
plot(mse ~ prop_yrs, xlab = "Proportion of extreme years observed",
     ylab = "Best model MSE", pch = 16, data = obs_clim_rng, col = "black")

dev.off()


# Models on climate/no -------------------------------------------------------------
mod_sample_size <- glm(model_climate ~ rep_n, family = "binomial", data = mod_climate)
mod_prop_rang   <- glm(model_climate ~ prop_rang, family = "binomial", data = mod_climate)
mod_prop_yrs    <- glm(model_climate ~ prop_yrs, family = "binomial", data = mod_climate)

mod_ecoregion   <- glm(model_climate ~ Ecoregion, family = "binomial", data = mod_climate)
mod_dicot_mono  <- glm(model_climate ~ DicotMonoc, family = "binomial", data = mod_climate)
mod_class       <- glm(model_climate ~ Class, family = "binomial", data = mod_climate)

# analyze results
res_summary_mod <- lapply( list(mod_sample_size, mod_prop_rang, mod_prop_yrs,
                            mod_ecoregion, mod_dicot_mono, mod_class) , summary) %>%
                      setNames(c("sample_size", "prop_rang", "prop_yrs",
                               "ecoregion", "dicot_mono","class"))
