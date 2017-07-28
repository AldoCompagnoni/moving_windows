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
pred_acc    <- Reduce(function(...) merge(...), pred_acc_l) %>%
                  rename( rep_n = rep_n_mse) %>%
                  inner_join( categ )

# calculate number of categories for each "best model"
best_mod_by_categ <- function(best_mod = NULL, category){
  
  if( is.null(best_mod) ){
    tmp <- pred_acc %>%
      group_by_( category ) %>%
      summarise( rep = n() ) %>%
      as.data.frame %>% t
  } else{
    tmp <- pred_acc %>%
      subset( model_mse == best_mod ) %>%
      group_by_( category ) %>%
      summarise( rep = n() ) %>%
      as.data.frame %>% t
  }
  df <- matrix(tmp[2,], nrow=1, ncol=length(tmp[2,]), byrow=T) %>%
    as.data.frame %>%
    setNames( tmp[1,] ) %>% as.matrix()
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
plot(mse ~ rep_n, pch = 16, data = pred_acc, 
     xlab = "Number of replicates", ylab = "Mean squared error")

dev.off()


# replication of best models
tiff(paste0("results/plots/best_mod_replication.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2), mar = c(3.5,3.5,1,0.5), mgp = c(2,0.7,0) ,
    cex.lab = 1.2)

hist(pred_acc$rep_n, ylim = c(0,8), xlim = c(0,70), xlab = "replication", main = "All data sets")
hist(subset(pred_acc, model_mse == "expp")$rep_n, ylim = c(0,8), xlim = c(0,70), xlab = "replication", main = "Power exp")
hist(subset(pred_acc, model_mse == "ctrl2")$rep_n, ylim = c(0,8), xlim = c(0,70), xlab = "replication", main = "24 month mean")
hist(subset(pred_acc, model_mse == "ctrl1")$rep_n, ylim = c(0,8), xlim = c(0,70), xlab = "replication", main = "NULL")

dev.off()


# best model by ecoregion
tiff(paste0("results/plots/best_mod_by_ecoregion.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2), mar = c(3,3.5,1.5,0.5), mgp = c(2,0.7,0),
    cex.lab = 1.2)

barplot( best_mod_by_categ(NULL,"Ecoregion"), main = "Representation across species", ylim = c(0,14))
barplot( best_mod_by_categ("expp","Ecoregion"), main = "Power Exponential", ylim = c(0,14) )
barplot( best_mod_by_categ("ctrl2","Ecoregion"), main = "24 Months", ylim = c(0,14) )
barplot( best_mod_by_categ("ctrl1","Ecoregion"), main = "NULL", ylim = c(0,14) )

dev.off()

# # best model by DicotMonoc
# tiff(paste0("results/plots/best_mod_by_DicotMonoc.tiff"),
#      unit="in", width=6.3, height=6.3, res=600,compression="lzw")
# 
# par(mfrow=c(2,2), mar = c(3,3.5,1.5,0.5), mgp = c(2,0.7,0),
#     cex.lab = 1.2)
# 
# barplot( best_mod_by_categ(NULL,"DicotMonoc"), main = "Representation across species", ylim = c(0,29))
# barplot( best_mod_by_categ("expp","DicotMonoc"), main = "Power Exponential", ylim = c(0,29) )
# barplot( best_mod_by_categ("ctrl2","DicotMonoc"), main = "24 Months", ylim = c(0,29) )
# barplot( best_mod_by_categ("ctrl1","DicotMonoc"), main = "NULL", ylim = c(0,29) )
# 
# dev.off()
# 
# # best model by Class
# tiff(paste0("results/plots/best_mod_by_Class.tiff"),
#      unit="in", width=6.3, height=6.3, res=600,compression="lzw")
# 
# par(mfrow=c(2,2), mar = c(3,3.5,1.5,0.5), mgp = c(2,0.7,0),
#     cex.lab = 1.2)
# 
# barplot( best_mod_by_categ(NULL,"Class"), main = "Representation across species", ylim = c(0,29))
# barplot( best_mod_by_categ("expp","Class"), main = "Power Exponential", ylim = c(0,29) )
# barplot( best_mod_by_categ("ctrl2","Class"), main = "24 Months", ylim = c(0,29) )
# barplot( best_mod_by_categ("ctrl1","Class"), main = "NULL", ylim = c(0,29) )
# 
# dev.off()





# best model by climate sampled --------------------------------------------------------------





# climate info ------------------------------------------------------------------------
spp_list      <- lambdas$SpeciesAuthor %>% unique %>% sort
spp_list      <- spp_list[-12] # what happened with Daphne_rodriguezii?!?

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



# best model by observed climate range
obs_clim_rng <- merge(pred_acc, clim_rng_m)

tiff(paste0("results/plots/best_mod_by_climate_sampled.tiff"),
     unit="in", width=3.5, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,1), mar = c(2,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1.2)

boxplot(prop_rang ~ model_mse, data = obs_clim_rng,
        names = c("Null", "24.mon", "Expp"), ylab = "Proportion of climate observed",
        ylab = "Proportion")
boxplot(prop_yrs ~ model_mse, data = obs_clim_rng,
        names = c("Null", "24.mon", "Expp"), ylab = "Proportion of extreme years",
        ylab = "Proportion")

dev.off()


par(mfrow=c(2,1), mar = c(3.5,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1.2)

plot(mse ~ prop_rang,xlab = "Proportion of climate observed", 
     pch = 16, data = obs_clim_rng)
plot(mse ~ prop_yrs, xlab = "Proportion of extreme years observed", 
     pch = 16, data = obs_clim_rng)


subset(obs_clim_rng, prop_rang == 1)
subset(obs_clim_rng, prop_yrs  == 0)
