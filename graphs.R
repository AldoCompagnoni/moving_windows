# first try to automate graphs 
setwd("C:/cloud/MEGA/projects/sAPROPOS/")
source("~/moving_windows/format_data.R")
library(dplyr)
library(tidyr)
options(stringsAsFactors = F )


# summary plots ---------------------------------------------------------------------
sum_files <- list.files("results")[grep("mod_summaries_", list.files("results") )]
crx_files <- list.files("results")[grep("crossval_", list.files("results") )]
mod_summ  <- lapply(sum_files, function(x) read.csv(paste0("results/",x)) ) %>%
                setNames( gsub("mod_summaries_", "", sum_files ) ) %>%
                setNames( gsub(".csv", "", names(.) ) )
replic    <- lapply(crx_files, function(x) read.csv(paste0("results/",x)) %>% nrow ) %>%
                setNames( gsub("crossval_pred_diag_", "", crx_files ) ) %>%
                setNames( gsub(".csv", "", names(.) ) ) %>%
                unlist %>% t %>% t %>% as.data.frame %>% 
                tibble::rownames_to_column(var = "species") %>%
                rename( rep_n = V1)

# best models
best_mods <- lapply(mod_summ, function(x) x[1,"model"]) %>%
                unlist %>%
                table
                
# predictive accuracy
pred_acc  <- lapply(mod_summ, function(x) x[1,c("mse")]) %>%
                unlist %>%
                t %>% t %>% 
                as.data.frame %>%
                tibble::rownames_to_column(var = "species") %>%
                rename( mse = V1 ) %>%
                inner_join( replic )


# the graph
tiff(paste0("results/plots/summaries.tiff"),
     unit="in", width=6.3,height=3.5,res=600,compression="lzw")

par(mfrow=c(1,2), mar = c(3.5,3.2,0.5,0.5), mgp = c(2,0.7,0) ,
    cex.lab = 1.2)
barplot(best_mods, ylab = "Best model", cex.names = 1.2)
plot(mse ~ rep_n, pch = 16, 
     data = pred_acc, 
     xlab = "Number of replicates", ylab = "Mean squared error")

dev.off()


# single species graphs --------------------------------------------------------------------
post_files  <- list.files("results")[grep("posterior_", list.files("results") )]

# create species list based on what's available!
spp_list    <- intersect(gsub(".csv", "", gsub("mod_summaries_", "", sum_files) ),
                         gsub(".csv", "", gsub("posterior_", "", post_files) )
                         )
  
# loop through (available) species
for(ii in 18:length(spp_list) ){
  
spp_name<- spp_list[ii] # test run w/ spp number 1
lam     <- read.csv("C:/cloud/Dropbox/sAPROPOS project/DemogData/lambdas_6tr.csv") 
m_info  <- read.csv("C:/cloud/MEGA/Projects/sApropos/MatrixEndMonth_information.csv")
clim    <- read.csv("C:/cloud/Dropbox/sAPROPOS project/DemogData/precip_fc_demo.csv") #%>%
              #mutate( ppt = ppt / 30)
spp     <- clim$species %>% unique

# add monthg info to lambda information
month_add <- m_info %>%
              mutate(SpeciesAuthor = trimws(SpeciesAuthor) ) %>%
              select(SpeciesAuthor, MatrixEndMonth)
lam_add   <- subset(lam, SpeciesAuthor %in% m_info$SpeciesAuthor) %>%  
              select(-MatrixEndMonth) %>%
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
post_df <- read.csv(paste0("C:/cloud/MEGA/Projects/sApropos/results/posterior_",
                    spp_name,".csv") ) %>%
              dplyr::select(-grep("log_lik", names(.)) ) %>%
              subset( model == "gaus" )
mean_df <- read.csv(paste0("C:/cloud/MEGA/Projects/sApropos/results/mod_summaries_",
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
                  select(-2) %>%
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
