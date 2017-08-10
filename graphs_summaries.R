# first try to automate graphs 
setwd("C:/cloud/MEGA/projects/sAPROPOS/")
source("C:/CODE/moving_windows/format_data.R")
library(dplyr)
library(tidyr)
library(testthat)
options(stringsAsFactors = F )

crossval_type <- "loyo"
clim_var      <- "airt"
m_back        <- 24


# summarize moving windows results by climate variable
summ_by_clim <- function(clim_var){
  
  # read lambda/clim data ---------------------------------------------------------------------------------
  lam     <- read.csv("lambdas_6tr.csv", stringsAsFactors = F)
  m_info  <- read.csv("MatrixEndMonth_information.csv", stringsAsFactors = F)
  clim    <- data.table::fread(paste0(clim_var,"_fc_demo.csv"), stringsAsFactors = F)
  
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
  mod_splin <- read.csv( paste0("results/splines/loyo/summaries/spline_",clim_var,"24_summaries.csv") ) %>%
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
    clim_separate <- clim_list(spp_name, clim, spp_lambdas)
    
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
                                     color    = c("blue", "orange", "red") ) )


# summary plots ----------------------------------------------------------------------------------

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

par(mfrow=c(2,1), mar = c(3.5,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1.2)

plot(mse ~ prop_rang,xlab = "Proportion of climate values observed",
     ylab = "Best model MSE", pch = 16, data = mw_summ_df, 
     col = mw_summ_df$color)
legend("topleft", c("Air temperature", "Potential evapotranspiration", "Precipitation"), pch = 16,
       col = unique(mw_summ_df$color), bty = "n")
plot(mse ~ prop_yrs, xlab = "Proportion of extreme years observed",
     ylab = "Best model MSE", pch = 16, data = mw_summ_df, 
     col = mw_summ_df$color)

dev.off()


# barplots ----------------------------------------------------------------

# color code for all (hopefully!) barplots
col_barplot <- c("red","red","orange","orange","blue","blue")

# Best models by replication
tiff(paste0("results/moving_windows/",crossval_type,"/plots/best_mod_replication.tiff"),
     unit="in", width=6.3, height=3.15, res=600,compression="lzw")

par(mfrow=c(1,2), mar = c(2.5,3,0.1,0.5), mgp = c(1.7,0.6,0),
    cex.axis=0.5, cex.lab = 1)
boxplot(rep_n ~ model_climate + clim_var, data=mw_summ_df, 
        ylab = "N. of reps (site-by-year)", names.cex=0.5,
        names= rep(c("NULL","Clim."),3), col = col_barplot )
boxplot(rep_yr ~ model_climate + clim_var, data=mw_summ_df, 
        ylab = "N. of years sampled",
        names= rep(c("NULL","Clim."),3), col = col_barplot )
legend("topright", c("Air temp.", "PET", "precip"),
       fill = unique(col_barplot), bty = "n")

dev.off()


# best models by climate sampled
tiff(paste0("results/moving_windows/",crossval_type,"/plots/best_mod_by_climate_sampled.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,2), mar = c(2,3.5,0.5,0.2), mgp = c(2,0.7,0), cex.lab = 1.2)

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
        ylab = "St. Dev. of mean sampled climate from mean hist. climate" )
abline(h = 0, lty = 2)
legend("topleft", c("Air temperature", "Potential evapotranspiration", "precipitation"),
       fill = unique(col_barplot), bty = "n")

dev.off()



# # best model by ecoregion
# tiff(paste0("results/moving_windows/",crossval_type,"/plots/",clim_var,"/best_mod_by_ecoregion(",clim_var,").tiff"),
#      unit="in", width=6.3, height=6.3, res=600,compression="lzw")
# 
# par(mfrow=c(2,2), mar = c(3,3.5,1.5,0.5), mgp = c(2,0.7,0), cex.lab = 1.2)
# barplot( best_mod_by_categ(NULL,"Ecoregion"), main = "Representation across species", 
#          ylim = c(0,14), col = "grey")
# barplot( best_mod_by_categ(1,"Ecoregion"), main = "Climate", ylim = c(0,14), col = "grey")
# barplot( best_mod_by_categ(0,"Ecoregion"), main = "NULL", ylim = c(0,14), col = "grey")
# 
# dev.off()


# graphs best spline models with/without climate ------------------------------------------------------
read_spline_summ <- function(clim_var, m_back){
  return(
    read.csv( paste0("results/splines/loyo/summaries/spline_",clim_var,m_back,"_summaries.csv") ) %>%
      dplyr::select(species,dev0,dev1) %>%
      unique %>%
      mutate( model_climate = 0 ) %>%
      mutate( model_climate = replace( model_climate, which(dev0 > dev1), 1) ) %>%
      setNames( c("species", "dev0_spline", "dev1spline", "model_climate_spline") ) %>%
      as.data.frame()
  )
  
}

# plots
tiff(paste0("results/splines/loyo/best_spline_models.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(2,3), mar = c(2,3.5,1.5,0.2), mgp = c(2,0.7,0) )

m_back <- 12
clim_var <- "precip"
barplot( table(read_spline_summ(clim_var,m_back)$model_climate_spline),
         names = c("NULL", clim_var),
         main = paste0(clim_var,", ",m_back," months") )
clim_var <- "pet"
barplot( table(read_spline_summ(clim_var,m_back)$model_climate_spline),
         names = c("NULL", "climate"),
         main = paste0(clim_var,", ",m_back," months") )
clim_var <- "airt"
barplot( table(read_spline_summ("airt",m_back)$model_climate_spline),
         names = c("NULL", "climate"),
         main = paste0(clim_var,", ",m_back," months") )

m_back <- 24
clim_var <- "precip"
barplot( table(read_spline_summ(clim_var,m_back)$model_climate_spline),
         names = c("NULL", clim_var),
         main = paste0(clim_var,", ",m_back," months") )
clim_var <- "pet"
barplot( table(read_spline_summ(clim_var,m_back)$model_climate_spline),
         names = c("NULL", "climate"),
         main = paste0(clim_var,", ",m_back," months") )
clim_var <- "airt"
barplot( table(read_spline_summ("airt",m_back)$model_climate_spline),
         names = c("NULL", "climate"),
         main = paste0(clim_var,", ",m_back," months") )

dev.off()





# Models on climate/no -------------------------------------------------------------
model_climate_mods <- list(
  model_climate ~ rep_n,
  model_climate ~ rep_yr,
  model_climate ~ prop_rang,
  model_climate ~ prop_yrs,
  model_climate ~ prop_var,
  model_climate ~ prop_var_r,
  model_climate ~ mean_dev,
  model_climate ~ mean_clim,
  model_climate ~ Ecoregion,
  model_climate ~ DicotMonoc,
  model_climate ~ Class
)

# fit models
models <- lapply(model_climate_mods, function(x) glm(x, family = "binomial", data = obs_clim_rng)) %>%
                    setNames(c("sample_size", "sampled_years", "prop_rang", "prop_yrs",
                               "mod_prop_var", "mod_prop_var_r", "mod_prop_mean","mod_mean_clim",
                               "ecoregion", "dicot_mono","class"))

# summarise model results
res_summary_mod <- lapply(models, summary)


# compare with splines --------------------------------------------------
compare_df <- merge(obs_clim_rng, mod_splin) %>%
                dplyr::select(species, model_climate, model_climate_spline)

mod <- glm(model_climate_spline ~ model_climate, family= "binomial", data = compare_df)
summary(mod)
