# produce summary graphs
setwd("C:/cloud/Dropbox/sAPROPOS/")
source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(magrittr)
library(testthat)
options(stringsAsFactors = F )

today_date    <- gsub("-","_",Sys.time() %>% as.Date %>% as.character)
m_back        <- 36
interval      <- NULL
pre_chelsa    <- NULL # '_pre_chelsa'

# pipe-able Reduce_rbind and grep functions
rbind_l     <- function(df_list){ Reduce(function(...) rbind(...), df_list) }


# summarize moving windows results by climate variable
mod_perform <- function(ii){
  
  clim_var <- input_df$clim_var[ii]
  response <- input_df$response[ii]
  interval <- input_df$interval[ii]
  resp_clim<- paste0("_",response,"_",clim_var)
  
  # read lambda/clim data ---------------------------------------------------------------------------------
  lam     <- read.csv("all_demog_6tr.csv", stringsAsFactors = F) #%>%
                #subset( SpeciesAuthor != "Purshia_subintegra" )

  # summary info ----------------------------------------------------------------------------
  res_folder<- paste0("results/moving_windows/",response,"/",clim_var,pre_chelsa,interval) 
  sum_files <- list.files(res_folder)[grep("mod_summaries_", list.files(res_folder) )] %>% 
                  stringr::str_subset(resp_clim)
  mod_summ  <- lapply(sum_files, function(x) read.csv(paste0(res_folder,"/",x)) ) %>%
                  setNames( gsub("mod_summaries_", "", sum_files ) ) %>%
                  setNames( gsub(paste0(resp_clim,".csv"), "", names(.) ) )
  all_sums  <- Map(function(x,y) tibble::add_column(x, species = y, .before = 1), 
                     mod_summ, names(mod_summ) ) %>% 
                  lapply(function(x) dplyr::select(x,species, model, waic, looic, mse, elpd)) %>% 
                  rbind_l  
          
  # species
  spp           <- names(mod_summ)
 
  # create information on total replication
  create_rep_n <- function(x, lam, response){
    
    data.frame( species = x,
                rep_n   = format_species(x, lam, response) %>% rbind_l %>% nrow,
                stringsAsFactors = F)
    
  }
  rep_n_df <- lapply(spp, create_rep_n, lam, response) %>% rbind_l
  
  # spit out 
  left_join(all_sums, rep_n_df) %>% 
    mutate(mean_elpd = elpd / rep_n,
           clim_var = clim_var,
           response = response)

}

# all models results
input_df    <- expand.grid( clim_var = c("precip","airt"),
                            response = c("surv","grow","fec","log_lambda"),
                            interval = "",
                            stringsAsFactors = F)

# ALL model information
mod_perf_df  <- lapply(1:nrow(input_df), mod_perform) %>%
                  rbind_l %>% 
                  arrange(species, mse)


# Matplot of relative (%) distances to best model
relative_bp <- function(resp, clim_v){
  
  perf_by_spp <- mod_perf_df %>% 
                    subset( response == resp & clim_var == clim_v ) %>% 
                    dplyr::select(species, model, mean_elpd) %>% 
                    mutate( mean_elpd = replace(mean_elpd, mean_elpd == -Inf, NA)) %>% 
                    spread(model, mean_elpd)
  perf_mat    <- dplyr::select(perf_by_spp, -species) %>% t
  
  # get minimum of maximum values
  get_min     <- function(ii) which(perf_mat[,ii] == min(perf_mat[,ii],na.rm=T))
  get_max     <- function(ii) which(perf_mat[,ii] == max(perf_mat[,ii],na.rm=T))
  
  # species sequence
  spp_seq     <- 1:ncol(perf_mat)
  
  # matrix of benchmark ids
  max_ids   <- sapply(spp_seq ,get_max) %>% 
                    cbind( spp_seq )
  min_ids   <- sapply(spp_seq, get_min) %>% 
                    cbind( spp_seq )
  
  # benchmark values
  max_val   <- perf_mat[max_ids]
  min_val   <- perf_mat[min_ids]
  scale_val <- max_val - min_val
  
  # differences matrix 
  diff_mat    <- sweep(perf_mat, 2,  max_val, FUN='-') %>% 
                    sweep(2, scale_val, FUN='/')*-1 
    
  # labels for models to plot on x-axis
  mod_labs <- c('ctrl1','ctrl2','yr1','yr2','yr3','yr_bet',
                'yr_wgt','gaus','expp','gev',
                'expp_n','gev_n','simpl_n')
  
  # differences data frame
  long_df <- diff_mat %>% 
                as.data.frame %>% 
                setNames( perf_by_spp$species ) %>% 
                tibble::add_column(model = row.names(.), .before=1) %>% 
                gather(species,mean_elpd,Actaea_spicata:Thelesperma_megapotamicum) %>% 
                mutate( model = factor(model, levels = mod_labs) )
            
  

  boxplot( mean_elpd ~ model, data = long_df,
         xaxt = "n",  xlab = "", ylab = "",
         main = paste0( resp, '~', clim_v) )
  
  # x axis with ticks but without labels
  axis(1, at=1:13, labels = FALSE)
    
  # Plot x labs at default x position
  text(x =  seq_along(mod_labs), 
       y = par("usr")[3] - 0.05, srt = 45, adj = 1,
       labels = mod_labs, xpd = TRUE)
      
}


tiff(paste0("results/moving_windows/prop_dist_to_best_model.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(4,2), mar= c(3, 1.5, 1, 0.2) + 0.1, oma = c(0,3.5,0,0),
    mgp = c(2,0.5,0))
relative_bp('surv','airt')
relative_bp('surv','precip')

relative_bp('grow','airt')
mtext("Proportion of distance from best model (mean LPPD)", 2,line=2.5, at = -0.3,
      cex = 1.3, xpd = NA)
relative_bp('grow','precip')

relative_bp('fec','airt')
relative_bp('fec','precip')

relative_bp('log_lambda','airt')
relative_bp('log_lambda','precip')

dev.off()



tiff(paste0("results/moving_windows/mean_lppd_by_resp_climvar.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

mod_perf_df <- mod_perf_df %>%
                  mutate( mean_elpd = replace(mean_elpd, mean_elpd == -Inf, NA) )

par(mfrow=c(1,1), mar = c(3,3,0.2,0.2), oma = c(0,0,0,0),
    mgp=c(1.5,0.5,0))
boxplot(mean_elpd ~  clim_var + response, data = mod_perf_df,
        ylab = "Mean LPPD", xaxt='n',cex.lab=1.3)
axis(1,at=1:8,labels=F)
text(1:8,y=-23.8,c('fec','fec','grow','grow','lambda','lambda','surv','surv'), xpd=NA)
text(1:8,y=-24.8,c('(airt)','(prec)','(airt)','(prec)','(airt)','(prec)','(airt)','(prec)'), xpd=NA)

abline(h=0,lty=2)

dev.off()


# absolute LPPD
abs_lppd <- function(resp, clim_v){
  
  # absolute measures
  boxp_df <- mod_perf_df %>% 
                subset( response == resp & clim_var == clim_v) %>% 
                mutate( mean_elpd = replace(mean_elpd,
                                            mean_elpd == -Inf,
                                            NA))
  
  if( grepl('grow|surv', resp) ){
  
    boxplot(mean_elpd ~ species, data=boxp_df,
            boxlwd = 0.5, medlwd = 1, outlwd = 0.5, whisklwd = 0.5,
            outcex = 0.5,
            xaxt = 'n', ylim = c(-3,6), cex.main = 1.5,
            main = paste0(resp, ' ~ ', clim_v) )
    # x axis with ticks but without labels
    axis(1, labels = F, tick=F)
    abline(h=0, lty = 2, lwd = 0.5)  
    abline(v=5, lty = 2, lwd = 0.5)  
    abline(v=15, lty = 2, lwd = 0.5)
    abline(v=25, lty = 2, lwd = 0.5)
  
  }else{
    
    boxplot(mean_elpd ~ species, data=boxp_df,
            boxlwd = 0.5, medlwd = 1, outlwd = 0.5, whisklwd = 0.5,
            outcex = 0.5,
            xaxt = 'n', ylim = c(-5,10), cex.main = 1.5,
            main = paste0(resp, ' ~ ', clim_v) )
    # x axis with ticks but without labels
    axis(1, labels = F, tick=F)
    abline(h=0, lty = 2, lwd = 0.5)  
    abline(v=5, lty = 2, lwd = 0.5)  
    abline(v=15, lty = 2, lwd = 0.5)
    abline(v=25, lty = 2, lwd = 0.5)

  }
  
}

# plot species labels
plot_spp_lab <- function(x){

  spp_labs <- unique(x$species) %>% 
                stringr::str_replace('var._gnaphalifolium_2','...')
  
  # Plot x labs at default x position
  text(x =  seq_along(spp_labs), 
       y = -6, srt = 90, adj = 1,
       labels = spp_labs, xpd = NA, cex = 0.7)
}



tiff(paste0("results/moving_windows/mean_lppd_by_spp.tiff"),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

par(mfrow=c(3,2), mar = c(0,2,2,0.2), oma = c(8,2,0,0),
      mgp=c(1.7,0.5,0))

abs_lppd('surv', 'airt')
abs_lppd('surv', 'precip')

abs_lppd('grow', 'airt')
mtext("Mean LPPD", 2,line=2, at = 3,
      cex = 1.3, xpd = NA)
abs_lppd('grow', 'airt')

abs_lppd('fec', 'airt')
plot_spp_lab(mod_perf_df)
abs_lppd('fec', 'airt')
plot_spp_lab(mod_perf_df)

dev.off()




# image_mat <- -(mod_sel_tab('grow','airt')[,-1] %>% as.matrix)
# image_mat <- -(mod_sel_tab('grow','precip')[,-1] %>% as.matrix)
# image_mat <- -(mod_sel_tab('surv','airt')[,-1] %>% as.matrix)
# image_mat <- -(mod_sel_tab('surv','precip')[,-1] %>% as.matrix)
# 
# image_mat <- -(mod_sel_tab('fec','airt')[,-1] %>% as.matrix)
# image_mat <- -(mod_sel_tab('fec','precip')[,-1] %>% as.matrix)
# 
# # image
# print(lattice::levelplot(perf_mat, col.regions=heat.colors))
# 
# 
# # image
# image(mod_sel_tab('grow','precip')[,-1] %>% as.matrix)
# 
# 
# 
# 
# # elpd heat map
# elpd_heat_map <- function(resp, clim_v){
#   
#   perf_by_spp <- mod_perf_df %>% 
#                     subset( response == resp & clim_var == clim_v ) %>% 
#                                   dplyr::select(species, model, mean_elpd) %>% 
#                                   spread(model, mean_elpd)
#   perf_mat    <- dplyr::select(perf_by_spp, -species) %>% t
#   perf_mat[perf_mat == -Inf] <- NA
#   
#   print(lattice::levelplot(t(perf_mat), col.regions=heat.colors))
# 
# }
# 
# 
# # ggplot heatmap --------------------------------------------------
# 
# # elpd heat map
# elpd_heat_map <- function(resp, clim_v){
#   
#   perf_by_spp <- mod_perf_df %>% 
#                   subset( response == resp & clim_var == clim_v ) %>% 
#                                 dplyr::select(model, species,  mean_elpd) %>%
#                   arrange( species, model) %>% 
#                   mutate( mean_elpd = replace(mean_elpd, mean_elpd == -Inf, NA) ) %>% 
#                   mutate( species = 
#                             replace(species, 
#                                     grepl('Eriogonum',species),
#                                     'Eriogonum_longifolium_var...') )
# 
#   
#   # plot it out
#   ggplot(perf_by_spp, aes(species, model )) +
#   geom_tile(aes(fill = mean_elpd), color = "white") +
#   scale_fill_gradient2(low = "white", high = "steelblue") +
#   ylab( paste0( 'Model: ',resp," ~ ",clim_v) ) +
#   xlab("Species ") +
#   theme(legend.title = element_text(size = 10),
#         legend.text = element_text(size = 12),
#         plot.title = element_text(size=16),
#         axis.title=element_text(size=14,face="bold"),
#         axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(fill = "Mean LPPD")
# 
# }
# 
# 
# # heat maps
# elpd_heat_map('surv','airt')
# elpd_heat_map('surv','precip')
# elpd_heat_map('grow','airt')
# elpd_heat_map('grow','precip')
# elpd_heat_map('fec','airt')
# elpd_heat_map('fec','precip')
# elpd_heat_map('log_lambda','airt')
# elpd_heat_map('log_lambda','precip')
# 
# 
# 
# 
# # matplot
# matplot_map <- function(resp, clim_v){
#   
#   perf_by_spp <- mod_perf_df %>% 
#                   subset( response == resp & clim_var == clim_v ) %>% 
#                                 dplyr::select(model, species,  mean_elpd) %>%
#                   arrange( species, model) %>% 
#                   mutate( mean_elpd = replace(mean_elpd, mean_elpd == -Inf, NA) ) %>% 
#                   mutate( species = 
#                             replace(species, 
#                                     grepl('Eriogonum',species),
#                                     'Eriogonum_longifolium_var...') ) %>% 
#                   spread(species,mean_elpd)
# 
#   matplot(perf_by_spp[,-1], type = 'l',
#           xlab = "model")
#   
# }
# 
# matplot_map('surv','airt')
# matplot_map('surv','precip')
# 
