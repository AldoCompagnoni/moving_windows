# first try to automate graphs 
setwd("C:/cloud/MEGA/projects/sAPROPOS/")
source("~/moving_windows/format_data.R")
library(dplyr)
library(tidyr)
options(stringsAsFactors = F )


# read species information ---------------------------------------------------------------------------------
lam     <- read.csv("C:/cloud/Dropbox/sAPROPOS project/DemogData/lambdas_6tr.csv", stringsAsFactors = F) 

# summaries organism type/biome
by_spp_rep <- function(var){
  
  lam %>%
    select( c("SpeciesAuthor", var) ) %>%
    distinct %>%
    group_by_(var) %>%
    summarise(spp_n = n())
  
}
# replication of categories
categ_rep <- lapply(c("OrganismType","Ecoregion","Order", "Class", "DicotMonoc", "Altitude"), by_spp_rep)

# I would use following categories       
categ     <- lam %>%
                select( c("SpeciesAuthor","Ecoregion", "DicotMonoc", "Class") ) %>%
                rename( species = SpeciesAuthor) %>%
                unique


# summary info ---------------------------------------------------------------------
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
  
  lapply(mod_summ, function(x) arrange_(x, mof)) %>%
    lapply(function(x) x[1,c(mof)]) %>%
    unlist %>% t %>% t %>% as.data.frame %>%
    tibble::rownames_to_column(var = "species") %>%
    inner_join( replic ) %>%
    setNames( c("species", mof, paste0("n_rep_",mof)) )
  
}

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


# summary plots ---------------------------------------------------------------------
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
     unit="in", width=6.3, height=3.15, res=600,compression="lzw")

par(mfrow=c(1,2), mar = c(3.5,3.2,0.5,0.5), mgp = c(2,0.7,0) ,
    cex.lab = 1.2)
plot(mse ~ rep_n, pch = 16, data = pred_acc, 
     xlab = "Number of replicates", ylab = "Mean squared error")
plot(looic ~ rep_n, pch = 16, data = pred_acc, 
     xlab = "Number of replicates", ylab = "looic")

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

# best model by DicotMonoc
tiff(paste0("results/plots/best_mod_by_DicotMonoc.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2), mar = c(3,3.5,1.5,0.5), mgp = c(2,0.7,0),
    cex.lab = 1.2)

barplot( best_mod_by_categ(NULL,"DicotMonoc"), main = "Representation across species", ylim = c(0,29))
barplot( best_mod_by_categ("expp","DicotMonoc"), main = "Power Exponential", ylim = c(0,29) )
barplot( best_mod_by_categ("ctrl2","DicotMonoc"), main = "24 Months", ylim = c(0,29) )
barplot( best_mod_by_categ("ctrl1","DicotMonoc"), main = "NULL", ylim = c(0,29) )

dev.off()

# best model by Class
tiff(paste0("results/plots/best_mod_by_Class.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2), mar = c(3,3.5,1.5,0.5), mgp = c(2,0.7,0),
    cex.lab = 1.2)

barplot( best_mod_by_categ(NULL,"Class"), main = "Representation across species", ylim = c(0,29))
barplot( best_mod_by_categ("expp","Class"), main = "Power Exponential", ylim = c(0,29) )
barplot( best_mod_by_categ("ctrl2","Class"), main = "24 Months", ylim = c(0,29) )
barplot( best_mod_by_categ("ctrl1","Class"), main = "NULL", ylim = c(0,29) )

dev.off()














# single species graphs --------------------------------------------------------------------
post_files  <- list.files(res_folder)[grep("posterior_", list.files(res_folder) )]

# create species list based on what's available!
spp_list    <- intersect(gsub(".csv", "", gsub("mod_summaries_", "", sum_files) ),
                         gsub(".csv", "", gsub("posterior_", "", post_files) )
                         )
  
# loop through (available) species
for(ii in 1:length(spp_list) ){
  
spp_name<- spp_list[ii] # test run w/ spp number 1
lam     <- read.csv("C:/cloud/Dropbox/sAPROPOS project/DemogData/lambdas_6tr.csv", stringsAsFactors = F) 
m_info  <- read.csv("C:/cloud/MEGA/Projects/sApropos/MatrixEndMonth_information.csv", stringsAsFactors = F)
clim    <- read.csv("C:/cloud/Dropbox/sAPROPOS project/DemogData/precip_fc_demo.csv", stringsAsFactors = F) #%>%
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


# plots ------------------------------------------------------------------------------------

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
