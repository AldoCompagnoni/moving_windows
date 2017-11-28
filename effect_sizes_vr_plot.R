# graphs of effect sizes
setwd("C:/cloud/MEGA/projects/sAPROPOS/")
source("C:/CODE/moving_windows/format_data.R")
library(tidyverse)
library(testthat)
options(stringsAsFactors = F )

m_back        <- 24
clim_var      <- c("airt","precip")
spp_cat       <- read.csv("sApropos_spp_raunkiaer.csv") %>% rename( species = SpeciesAuthor )
sens_df       <- read.csv("allCOMPADRE_COMADREOutput.csv") %>%
                    #subset( MatrixComposite == "Individual" &
                    subset( MatrixComposite == "Mean" &
                            Treatment == "Unmanipulated" ) %>%
                    rename( species = SpeciesAuthor,
                            surv    = Ssurv,
                            grow    = Sgrow,
                            fec     = Srep) %>%
                    dplyr::select( species, surv, grow, fec) %>%
                    group_by( species ) %>%
                    summarize_all( mean, na.rm=T ) %>%
                    subset( !(is.nan(surv) & is.nan(grow) & is.nan(fec)) ) %>%
                    gather(response, sensitivity, -species )

# create each files' path
path_create <- function(x,y){
  post_files <- Filter(function(x) grepl("posterior",x), x)
  paste0(y,post_files)
}

# calculate posterior means
post_means <- function(x){

  # extract species name
  spp_n     <- gsub( paste0("results/moving_windows/",response,"/[a-z]{3,6}/posterior_"),"",x)
  spp_name  <- gsub(".csv","",spp_n)

  # extract climate variable name
  # substitute
  # 1 everything before "/posterior_"
  # 2 everything after  "response/"
  clim_var  <- gsub(paste0("^(.*?)/",response,"/|/posterior_?.*"),"",x)

  posterior <- data.table::fread( paste0(x) ) %>%
    group_by( model ) %>%
    summarise_all( mean ) %>%
    ungroup %>%
    dplyr::select( -grep("log_lik_",names(.)) ) %>%
    as.data.frame %>%
    mutate(species = spp_name,
           clim_var = clim_var )

}

# single response variables -----------------------------------------------------------

# growth 
response      <- "grow"
wds           <- paste0("results/moving_windows/",response,"/",clim_var,"/") %>% as.list
raw_file_vec  <- lapply(wds, list.files)

# paths of files containing posteriors
post_paths_l  <- Map(path_create, raw_file_vec, wds)
post_paths    <- Reduce(function(...) c(...), post_paths_l) %>% as.list
means_l       <- lapply(post_paths, post_means)
means_g       <- Reduce(function(...) bind_rows(...), means_l) %>%
                    mutate( response = response)

# survival 
response      <- "surv"
wds           <- paste0("results/moving_windows/",response,"/",clim_var,"/") %>% as.list
raw_file_vec  <- lapply(wds, list.files)

# paths of files containing posteriors
post_paths_l  <- Map(path_create, raw_file_vec, wds)
post_paths    <- Reduce(function(...) c(...), post_paths_l) %>% as.list
means_l       <- lapply(post_paths, post_means)
means_s       <- Reduce(function(...) bind_rows(...), means_l) %>%
                    mutate( response = response)

# fecundity
response      <- "fec"
wds           <- paste0("results/moving_windows/",response,"/",clim_var,"/") %>% as.list
raw_file_vec  <- lapply(wds, list.files)

# paths of files containing posteriors
post_paths_l  <- Map(path_create, raw_file_vec, wds)
post_paths    <- Reduce(function(...) c(...), post_paths_l) %>% as.list
means_l       <- lapply(post_paths, post_means)
means_f       <- Reduce(function(...) bind_rows(...), means_l) %>%
                    mutate( response = response)

# rho 
response      <- "rho"
wds           <- paste0("results/moving_windows/",response,"/",clim_var,"/") %>% as.list
raw_file_vec  <- lapply(wds, list.files)

# paths of files containing posteriors
post_paths_l  <- Map(path_create, raw_file_vec, wds)
post_paths    <- Reduce(function(...) c(...), post_paths_l) %>% as.list
means_l       <- lapply(post_paths, post_means)
means_rho     <- Reduce(function(...) bind_rows(...), means_l) %>%
                    mutate( response = response)

# react_fsa (transient dynamics)
response      <- "react_fsa"
wds           <- paste0("results/moving_windows/",response,"/",clim_var,"/") %>% as.list
raw_file_vec  <- lapply(wds, list.files)

# paths of files containing posteriors
post_paths_l  <- Map(path_create, raw_file_vec, wds)
post_paths    <- Reduce(function(...) c(...), post_paths_l) %>% as.list
means_l       <- lapply(post_paths, post_means)
means_t       <- Reduce(function(...) bind_rows(...), means_l) %>%
                    mutate( response = response)

# react_fsa (transient dynamics)
response      <- "log_lambda"
wds           <- paste0("results/moving_windows/",response,"/",clim_var,"/") %>% as.list
raw_file_vec  <- lapply(wds, list.files)

# paths of files containing posteriors
post_paths_l  <- Map(path_create, raw_file_vec, wds)
post_paths    <- Reduce(function(...) c(...), post_paths_l) %>% as.list
means_l       <- lapply(post_paths, post_means)
means_ll      <- Reduce(function(...) bind_rows(...), means_l) %>%
                    mutate( response = response)


#  put it all together --------------------------------------------------------------------

# put it all together
means_l       <- list(means_s, means_g, means_f, means_rho, means_t, means_ll)
means_df_orig <- Reduce(function(...) rbind(...), means_l) %>%
                    mutate( response = as.factor(response) ) %>%
                    left_join( spp_cat ) %>%
                    mutate( col = "blue" ) %>%
                    mutate( col = replace(col, GrowthFormRaunkiaer == "Geophyte" | GrowthFormRaunkiaer == "Geophyte", "brown") ) %>%
                    mutate( col = replace(col, GrowthFormRaunkiaer == "Therophyte", "green") ) %>%
                    mutate(beta_abs = abs(beta) )

# REMOVE PROBLEMATIC SPECIES (need to model-check these one by one)
means_df_clean <- subset(means_df_orig, species != "Eryngium_alpinum" &
                                        species != "Psoralea_tenuiflora" & 
                                        species != "Purshia_subintegra" )

# add sensitivities and calculate sensitivity to climate
means_df       <- means_df_clean %>% 
                    left_join( sens_df ) %>%
                    mutate( clim_sens       = ( beta * exp(alpha) ) / ( (exp(alpha) +1)^2 ) ) %>%
                    mutate( clim_sens       = replace(clim_sens, 
                                                      response == "fec", 
                                                      beta[response == "fec"] * exp(alpha[response == "fec"])) ) %>%
                    mutate( clim_sens_abs = abs(clim_sens) ) %>%
                    mutate( clim_sens_scale = clim_sens * sensitivity ) %>%
                    mutate( clim_sens_scale_abs = abs(clim_sens_scale) )


# plot beta by response variable --------------------------------------------------------

# t_df        <- subset(means_df, response == "rho" | response == "react_fsa" )
# colors_df   <- dplyr::select(t_df, species, col) %>% unique %>% arrange()
# 
# # produce data frames
# t_beta_prec_mat <- subset(t_df, model == "ctrl2" & clim_var == "precip") %>%
#                       dplyr::select( species, response, beta_abs) %>%
#                       spread( species, beta_abs ) %>%
#                       .[c("response", colors_df$species)] %>%
#                       mutate( response = c(1,2) )
# t_beta_airt_mat <- subset(t_df, model == "ctrl2" & clim_var == "airt") %>%
#                       dplyr::select( species, response, beta_abs) %>%
#                       spread( species, beta_abs ) %>%
#                       .[c("response", colors_df$species)] %>%
#                       mutate( response = c(1,2) )
# 
# 
# # slopes, by vital rate and species
# tiff("results/moving_windows/slope_by_spp_transients.tiff", 
#      unit="in", width=3.15, height=6.3, res=400, compression="lzw")
# 
# par(mfrow=c(2,1), mar = c(2,2.5,1,1), mgp = c(1.5, 0.7,0) )
# 
# # Air temperature
# matplot(t_beta_airt_mat[,1], t_beta_airt_mat[,-1], 
#         type = "l", lty = 1, lwd = 2,
#         xlim = c(0.95, 2.05), ylim = c(0, 5), col = colors_df$col,
#         xaxt="n", xlab = "", ylab = expression(beta*" value"), main = "Air temperature")
# axis(1, at = c(1,2), labels = c("react_fsa", "rho") )
# legend("topleft", c("Hemicryptophyte + Nanophanerophyte","Geophyte + Chamaephyte","Therophyte"), 
#        lty = 1, lwd = 2, col = c("blue", "brown", "green" ), bty = "n", cex = 0.7)
# 
# # Precipitation
# matplot(t_beta_prec_mat[,1], t_beta_prec_mat[,-1], 
#         type = "l", lty = 1, lwd = 2,  
#         xlim = c(0.95, 2.05), ylim = c(0, 5), col = colors_df$col,
#         xaxt="n", xlab = "", ylab = expression(beta*" value"), main = "Precipitation")
# axis(1, at = c(1,2), labels = c("react_fsa", "rho") )
# 
# dev.off()



# climate_sensitivity by vital rate --------------------------------------------------------

vr_df         <- subset(means_df, response == "fec" | response == "grow" | response == "surv")
colors_df     <- dplyr::select(vr_df, species, col) %>% unique %>% arrange()

# produce data frames
vr_beta_prec_mat <- subset(vr_df, model == "ctrl2" & clim_var == "precip") %>%
                      dplyr::select( species, response, clim_sens_scale_abs) %>%
                      spread( species, clim_sens_scale_abs ) %>%
                      .[c("response", colors_df$species)] %>%
                      mutate( response = as.factor(as.character(response)) )

vr_beta_airt_mat <- subset(vr_df, model == "ctrl2" & clim_var == "airt") %>%
                      dplyr::select( species, response, clim_sens_scale_abs) %>%
                      spread( species, clim_sens_scale_abs ) %>%
                      .[c("response", colors_df$species)] %>%
                      mutate( response = as.factor(as.character(response)) )


# slopes, by vital rate and species
tiff("results/moving_windows/clim_sens_by_spp_vr.tiff",
     unit="in", width=6.3, height=3.5, res=400, compression="lzw")

par(mfrow=c(1,2), mar = c(2,2.5,1,1), mgp = c(1.5, 0.7,0) )

# # Air temperature
matplot(as.numeric(vr_beta_airt_mat[,1]), vr_beta_airt_mat[,-c(1,18)],
        type = "l", lty = 1, lwd = 2,
        xlim = c(0.9, 3.1), ylim = c(0, 1), col = colors_df$col[-17],
        xaxt="n", xlab = "", ylab = expression(italic(Sensitivity)*"  to Air Temperature"))
axis(1, at = c(1,2,3), labels=c("fecundity","growth", "survival") )
legend("topleft", c("Hemicryptophyte + Nanophanerophyte","Geophyte + Chamaephyte","Therophyte"),
       lty = 1, lwd = 2, col = c("blue", "brown", "green" ), bty = "n", cex = 0.7)

# # Precipitation
matplot(as.numeric(vr_beta_prec_mat[,1]), vr_beta_prec_mat[,-c(1,18)],
        type = "l", lty = 1, lwd = 2,
        xlim = c(0.9, 3.1), ylim = c(0, 1), col = colors_df$col[-17],
        xaxt="n", xlab = "", ylab = expression(italic(Sensitivity)*"  to Precipitation") )
axis(1, at = c(1,2,3), labels=c("fecundity","growth", "survival") )

dev.off()


tiff("results/moving_windows/sens_by_clim.tiff",
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

sens_plot <- data.frame( AirTemperature = stack(vr_beta_airt_mat[,-1])[,1],
                         Precipitation  = stack(vr_beta_prec_mat[,-1])[,1] ) %>%
                gather(clim_var, sens, AirTemperature, Precipitation)

par(mfrow=c(1,1), mar = c(2.5,3.5,0.1,0.1), mgp = c(2,0.7,0) )
boxplot(sens ~ clim_var, data = subset(sens_plot, sens < 12),
        names = c("ciao", "salve"),
        ylab = "Sensitivity", cex.lab = 2, xaxt = "n")
axis(1, at = c(1,2), c("Air temperature", "Precipitation"), cex.axis = 1.5)

dev.off()


# slopes, by vital rate 
tiff("results/moving_windows/slope_by_response.tiff", 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par(mfrow=c(2,1), mar = c(2,2.5,1,1), mgp = c(1.5, 0.5,0) )

boxplot(beta ~ response, subset(means_df, clim_var == "airt" & model == "ctrl2" & beta_abs < 4),
        main = "Air temperature", ylab = expression("Scaled "*beta*" value"),
        names = c("fec","grow",expression("log("*lambda*")"),"react_fsa","rho","surv") )
abline( h = 0, lwd=2 , lty = 2)
boxplot(beta ~ response, subset(means_df, clim_var == "precip" & model == "ctrl2"),
        main = "Precipitation", ylab = expression("Scaled "*beta*" value"),
        names = c("fec","grow",expression("log("*lambda*")"),"react_fsa","rho","surv") )
abline( h = 0, lwd=2 , lty = 2)

dev.off()


# slopes, by vital rate 
tiff("results/moving_windows/stzed_slope_by_response.tiff", 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par(mfrow=c(2,1), mar = c(2,2.5,1,1), mgp = c(1.5, 0.5,0) )

boxplot(clim_sens ~ response, subset(means_df, clim_var == "airt" & model == "ctrl2"),
        main = "Air temperature", ylab = expression("Scaled "*beta*" value") )
boxplot(abs(clim_sens) ~ response, subset(means_df, clim_var == "precip" & model == "ctrl2" & abs(clim_sens) < 5),
        main = "Precipitation", ylab = expression("Scaled "*beta*" value") )

dev.off()



# climate effect-by-sensitivity plots ------------------------------------------------------- 

# mean sensitivities by species
all_d <- read.csv("allCOMPADRE_COMADREOutput.csv") %>%
            subset( MatrixComposite == "Individual" &
                      Treatment == "Unmanipulated" & 
                      SpeciesAuthor %in% unique(means_df$species) ) %>%
            rename( species = SpeciesAuthor,
                    surv    = Ssurv,
                    grow    = Sgrow,
                    fec     = Srep) %>%
            dplyr::select( species, surv, grow, fec) %>%
            group_by( species ) %>%
            summarize_all( mean, na.rm=T ) %>%
            subset( !(is.nan(surv) & is.nan(grow) & is.nan(fec)) ) %>%
            gather(response, sensitivity, -species )


# climate effects + sensitivities 
# means_s_df <- means_df %>% left_join( all_d )
means_s_df <- means_df %>% left_join( sens_df )


beta_vs_sens <- means_s_df %>% 
                  subset( model == "ctrl2" & species != "Kosteletzkya_pentacarpos") %>%
                  subset(!(response == "log_lambda" | response == "rho" | response == "react_fsa") )
bvs_airt <- subset(beta_vs_sens, clim_var == "airt")
bvs_prec <- subset(beta_vs_sens, clim_var == "precip")


tiff("results/moving_windows/slope_by_vr_sens.tiff", 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par(mfrow = c(1,1) , mar = c(3,3,0.2,0.2), mgp = c(2,1,0) )
plot( clim_sens_abs ~ sensitivity, data=beta_vs_sens, pch = 16,
      ylab = expression("climate sensitivity"), col = as.factor(beta_vs_sens$response),
      xlab = expression("Sensitivity of "*lambda*" to vital rate") )
legend("topright", c("fecundity", "growth", "survival"), 
       col = c(1:3), pch = 16, bty = "n")

dev.off()


tiff("results/moving_windows/slope_by_vr_sens(stdz)_log.tiff", 
     unit="in", width=4, height=6.3, res=400, compression="lzw")

par(mfrow = c(2,1) , mar = c(3,3,1,0.2), mgp = c(2,1,0) )
plot( log(clim_sens_abs) ~ log(sensitivity), data= bvs_airt, pch = 16,
      ylab = expression("log[Sensitivity to temperature]"), col = as.factor(bvs_airt$response),
      xlab = expression("log[Sensitivity of "*lambda*" to vital rate]"),
      main = "Air temperature")
legend("bottomleft", c("fecundity", "growth", "survival"), 
       col = c(1:3), pch = 16, bty = "n")

plot( log(clim_sens_abs) ~ log(sensitivity), data= bvs_prec, pch = 16,
      ylab = expression("log[Sensitivity to precipitation]"), col = as.factor(bvs_prec$response),
      xlab = expression("log[Sensitivity of "*lambda*" to vital rate]"),
      main = "Precipitation")

dev.off()


# tiff("results/moving_windows/slope_by_vr_sens_log.tiff", 
#      unit="in", width=4, height=6.3, res=400, compression="lzw")
# 
# par(mfrow = c(2,1) , mar = c(3,3,1,0.2), mgp = c(2,1,0) )
# plot( log(beta_abs) ~ log(sensitivity), data= bvs_airt, pch = 16,
#       ylab = expression("log[climate effect ("*beta*")]"), col = as.factor(bvs_airt$response),
#       xlab = expression("log[Sensitivity of "*lambda*" to vital rate]"),
#       main = "Air temperature")
# legend("bottomleft", c("fecundity", "growth", "survival"), 
#        col = c(1:3), pch = 16, bty = "n")
# 
# plot( log(beta_abs) ~ log(sensitivity), data= bvs_prec, pch = 16,
#       ylab = expression("log[climate effect ("*beta*")]"), col = as.factor(bvs_prec$response),
#       xlab = expression("log[Sensitivity of "*lambda*" to vital rate]"),
#       main = "Precipitation")
# 
# dev.off()









# produce data frames -------------------------------------------------------------------------

# precipitation
tmp_prec_mat      <- subset(means_df, model == "ctrl2" & clim_var == "precip") %>%
  dplyr::select( species, response, beta ) %>%
  spread( species, beta ) %>% 
  t
vr_beta_prec_mat  <- tmp_prec_mat[-1,] %>%
  as.data.frame %>%
  apply(2, as.numeric) %>%
  as.data.frame %>%
  add_column( species = rownames(tmp_prec_mat)[-1] , .before = 1) %>%
  setNames( c("species", "fec_p", "grow_p", "surv_p") )

# temperature
tmp_airt_mat      <- subset(means_df, model == "ctrl2" & clim_var == "airt") %>%
  dplyr::select( species, response, beta ) %>%
  spread( species, beta ) %>% 
  t
vr_beta_airt_mat  <- tmp_airt_mat[-1,] %>%
  as.data.frame %>%
  apply(2, as.numeric) %>%
  as.data.frame %>%
  add_column( species = rownames(tmp_airt_mat)[-1] , .before = 1) %>%
  setNames( c("species", "fec_t", "grow_t", "surv_t") ) %>%
  left_join( colors_df )

# all betas
vr_betas          <- full_join(vr_beta_prec_mat, vr_beta_airt_mat)
write.csv(vr_betas, "sApropos_betas.csv", row.names=F)



# VR vs. VR (by temperature/precip)
tiff("results/moving_windows/vr_vs_by_climate.tiff", 
     unit="in", width=5.5, height=8.9, res=400, compression="lzw")

par( mfcol = c(3,2), mar = c(3,3,1,0.1), mgp = c(1.5, 0.5, 0), cex = 1.1)

# precipitation ---------------------------------------------

# growth vs. fecundity
plot(vr_betas$grow_p, vr_betas$fec_p,
     main = "Precipitation", pch = 16, col = vr_betas$col, 
     xlab = expression("Scaled "*beta*": Growth"),
     ylab = expression("Scaled "*beta*": Fecundity") )
cor_val <- cor( na.omit(vr_betas[,c("grow_p","fec_p")]) )[1,2]
legend("topright", legend = round(cor_val, 3), bty = "n", cex = 0.7)
p_val   <- summary( lm(vr_betas$fec_p ~ vr_betas$grow_p) )$coefficient[2,4]
legend("bottomright", legend = paste0("P: ", round(p_val, 3)), bty = "n", cex = 0.7)
legend("topleft", c("Hemicryptophyte","Geophyte","Therophyte"), 
       col = c("Blue","Brown","Green"), pch = 16, bty = "n", cex = 0.7)

# growth vs. survival
plot(vr_betas$grow_p, vr_betas$surv_p, 
     pch = 16, col = vr_betas$col,
     xlab = expression("Scaled "*beta*": Growth"),
     ylab = expression("Scaled "*beta*": Survival") )
cor_val <- cor( na.omit(vr_betas[,c("grow_p","surv_p")]) )[1,2]
legend("topright", legend = round(cor_val, 3), bty = "n", cex = 0.7)
p_val   <- summary( lm(vr_betas$surv_p ~ vr_betas$grow_p) )$coefficient[2,4]
legend("bottomright", legend = paste0("P: ", round(p_val, 3)), bty = "n", cex = 0.7)

# survival vs. fecundity
plot(vr_betas$surv_p, vr_betas$fec_p, 
     pch = 16, col = vr_betas$col,
     xlab = expression("Scaled "*beta*": Survival"),
     ylab = expression("Scaled "*beta*": Fecundity") )
cor_val <- cor( na.omit(vr_betas[,c("surv_p","fec_p")]) )[1,2]
legend("topright", legend = round(cor_val, 3), bty = "n", cex = 0.7)
p_val   <- summary( lm(vr_betas$fec_p ~ vr_betas$surv_p) )$coefficient[2,4]
legend("bottomright", legend = paste0("P: ", round(p_val, 3)), bty = "n", cex = 0.7)


# temperature ---------------------------------------------

# growth vs. fecundity
plot(vr_betas$grow_t, vr_betas$fec_t,
     main = "Temperature", pch = 16, col = vr_betas$col,
     xlab = expression("Scaled "*beta*": Growth"),
     ylab = expression("Scaled "*beta*": Fecundity") )
cor_val <- cor( na.omit(vr_betas[,c("grow_t","fec_t")]) )[1,2]
legend("topright", legend = round(cor_val, 3), bty = "n", cex = 0.7)
p_val   <- summary( lm(vr_betas$fec_t ~ vr_betas$grow_t) )$coefficient[2,4]
legend("bottomright", legend = paste0("P: ", round(p_val, 3)), bty = "n", cex = 0.7)

# growth vs. survival
plot(vr_betas$grow_t, vr_betas$surv_t, pch = 16,
     xlab = expression("Scaled "*beta*": Growth"), col = vr_betas$col,
     ylab = expression("Scaled "*beta*": Survival") )
cor_val <- cor( na.omit(vr_betas[,c("grow_t","surv_t")]) )[1,2]
legend("topright", legend = round(cor_val, 3), bty = "n", cex = 0.7)
p_val   <- summary( lm(vr_betas$surv_t ~ vr_betas$grow_t) )$coefficient[2,4]
legend("bottomright", legend = paste0("P: ", round(p_val, 3)), bty = "n", cex = 0.7)

# survival vs. fecundity
plot(vr_betas$surv_t, vr_betas$fec_t, pch = 16,
     xlab = expression("Scaled "*beta*": Survival"), col = vr_betas$col,
     ylab = expression("Scaled "*beta*": Fecundity") )
cor_val <- cor( na.omit(vr_betas[,c("surv_t","fec_t")]) )[1,2]
legend("topright", legend = round(cor_val, 3), bty = "n", cex = 0.7)
p_val   <- summary( lm(vr_betas$fec_t ~ vr_betas$surv_t) )$coefficient[2,4]
legend("bottomright", legend = paste0("P: ", round(p_val, 3)), bty = "n", cex = 0.7)

dev.off()


# Temperature vs. precip (by VR)
tiff("results/moving_windows/airt_v_prec_by_vr.tiff", 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par( mfrow = c(2,2), mar = c(3,3,1,0.1), mgp=c(1.7,0.5,0), oma = c(0,0,0,0.5) )

# fecundity
plot(vr_betas$fec_t, vr_betas$fec_p,
     xlab = expression("Scaled "*beta*": Temperature"),
     ylab = expression("Scaled "*beta*": Precipitation"), 
     main = "Fecundity", col = vr_betas$col, pch = 16 )

# growth
plot(vr_betas$grow_t, vr_betas$grow_p,
     xlab = expression("Scaled "*beta*": Temperature"),
     ylab = expression("Scaled "*beta*": Precipitation"), 
     main = "Growth",
     pch = 16, col = vr_betas$col )
legend("topleft", c("Hemicryptophyte","Geophyte","Therophyte"), 
       col = c("Blue","Brown","Green"), pch = 16, bty = "n", cex = 1)

# survival
plot(vr_betas$surv_t, vr_betas$surv_p,
     xlab = expression("Scaled "*beta*": Temperature"),
     ylab = expression("Scaled "*beta*": Precipitation"), 
     main = "Survival",
     pch = 16, col = vr_betas$col )

dev.off()


# Beta_prec by Beta_airt   
tiff("results/moving_windows/beta_tempPrec_line_vr.tiff", 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par( mfrow=c(1,1) )
ii = 1
plot(as.numeric(vr_betas[ii,2:4]),
     as.numeric(vr_betas[ii,5:7]),
     ylim = c(-4,4), xlim = c(-4,4), type = "b", pch = c(1:3), lwd = 2,
     xlab = expression(beta* " (Precipitation)"), 
     ylab = expression(beta* " (Temperature)"), col = vr_betas[ii,"col"] )

plot_linez <- function(xx){
  lines(as.numeric(vr_betas[xx,2:4]), type = "b", lwd = 2, pch = c(1:3),
        as.numeric(vr_betas[xx,5:7]), col = vr_betas[xx,"col"] )
}
lapply(2:nrow(vr_betas), plot_linez)
 
legend("topleft", c("fecundity", "growth", "survival"),
       pch = c(1:3), bty = "n")
legend("bottomleft", c("hemi+nano", "geo+chamae", "thero"),
       col = c("black", "grey", "brown"), lty = 1, lwd = 2, bty = "n")

dev.off()


# write Betas for Rob
subset(means_df, model == "ctrl2" & beta_abs < 5) %>%
  dplyr::select(species, clim_var, response, beta) %>%
  rename( SpeciesAuthor = species) %>%
  write.csv("results/moving_windows/betas_by_sppVr_sApropos.csv", row.names=F)
