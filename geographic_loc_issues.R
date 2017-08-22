setwd("C:/cloud/Dropbox/sAPROPOS project/")
library(tidyverse)

# read data -----------------------------------------------------------------------------------------
lam     <- read.csv("DemogData/lambdas_6tr.csv", stringsAsFactors = F) 
m_info  <- read.csv("C:/cloud/MEGA/Projects/sApropos/MatrixEndMonth_information.csv", stringsAsFactors = F)
clim    <- read.csv("DemogData/precip_fc_demo.csv",  stringsAsFactors = F) #%>%
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

mod_sum_l <- list()
# spp list - WITHOUT Daphne rodriguezii
spp_list  <- lambdas$SpeciesAuthor %>% unique %>% .[-3] # remove Daphne rodriguezii



# CLIMATE across sites (yearly means, yearly and monthly correlations) ---------------------------------------------

# identify species with multiple sites 
spp_mpl_sites <- function(ii){
  
  spp_name      <- spp_list[ii] # test run w/ spp number 1
  m_back        <- 24     # months back
  expp_beta     <- 20
  
  # lambda data
  spp_lambdas   <- format_species(spp_name, lambdas)
  
  data.frame( species = spp_name,
              site_rep = length(spp_lambdas) )  
  
}

# mean yearly climate: spp-by-site 
spp_by_site_climate <- function(ii){
  
  spp_name      <- spp_list[ii] # test run w/ spp number 1
  m_back        <- 24     # months back
  expp_beta     <- 20
  
  # lambda data
  spp_lambdas   <- format_species(spp_name, lambdas)
  
  # climate daily data
  clim_separate <- clim_list(spp_name, clim, spp_lambdas)
  
  # get yearly sums
  clim_sep_mean <- lapply(clim_separate, function(x)
                                            group_by(x, species, year ) %>%
                                            summarise( ppt = sum(ppt) ) %>%
                                            ungroup
                          ) 
  
  # get yearly mean ppt
  clim_means_l  <- Map(function(x,y) data.frame( species= spp_name,
                                                 site   = y,
                                                 mean_p = mean(x$ppt), 
                                                 sd_p   = sd(x$ppt) ),
                                      clim_separate, 
                                      names(clim_separate) 
                       )
  
  if( length(clim_means_l) > 1 ) {
    
    clim_means <- Reduce(function(...) rbind(...), clim_means_l)
     
  } else { 
    
    clim_means <- clim_means_l[[1]]
    
  }
  
  for(i in 1:2) clim_means[,i] <- as.character(clim_means[,i]) 
  clim_means
  
}


# between site differences in climate 
bw_site_clim <- function(ii){
  
  spp_name    <- mpl_site_spp[ii]
  
  sites_clim  <- subset(site_clim_df, species %in% spp_name) 
  cl_st_diff  <- data.frame( species = spp_name,
                             site_clim_diff = diff(sites_clim$mean_p) %>% as.matrix %>% as.vector
  )
  
}

# correlation of YEARLY climates between sites 
site_yr_clim_corr <- function(spp_name){
  
  #spp_name      <- spp_list[ii] # test run w/ spp number 1
  m_back        <- 24     # months back
  expp_beta     <- 20
  
  # lambda data
  spp_lambdas   <- format_species(spp_name, lambdas)
  
  # climate daily data
  clim_separate <- clim_list(spp_name, clim, spp_lambdas)
  
  # get yearly sums
  clim_sep_mean <- lapply(clim_separate, function(x)
    group_by(x, species, year ) %>%
      summarise( ppt = sum(ppt) ) %>%
      ungroup
  ) 
  
  # get correlation of climates across sites 
  site_yr_clim_l    <- lapply(clim_sep_mean, function(x) x$ppt )
  site_yr_clim_mat  <- do.call(cbind, site_yr_clim_l)
  cor_mat           <- cor(site_yr_clim_mat) 
  
  cor_mat[ upper.tri(cor_mat) ]
  
}

# correlation of MONTHLY climates between sites 
site_mon_clim_corr <- function(spp_name){
  
  #spp_name      <- spp_list[ii] # test run w/ spp number 1
  m_back        <- 24     # months back
  expp_beta     <- 20
  
  # lambda data
  spp_lambdas   <- format_species(spp_name, lambdas)
  
  # climate daily data
  clim_separate <- clim_list(spp_name, clim, spp_lambdas)
  
  
  # extractly monthly climate population-level climate; put it in "long" form
  clim_month <- function(clim_x, clim_var = "precip"){ #, pops
    
    # format day one
    day_one   <- as.Date( paste0("1/1/", first(clim_x$year) ), 
                          format="%d/%m/%Y")
    
    # climate data
    clim_d    <- as.Date(1:nrow(clim_x), day_one-1) %>%
      as.character %>%
      as.data.frame(stringsAsFactors=F) %>%
      separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
      bind_cols(clim_x) %>%
      dplyr::select(-year,-day) %>%
      setNames( c("year", "month", "day", "species", "ppt") )
    
    # monthly climates
    # if climate var relate to precipitation, monthly SUMS 
    if( clim_var == "precip" | clim_var == "pet"){
      clim_m   <- clim_d %>%
        group_by(year, month) %>%
        summarise( ppt = sum(ppt, na.rm=T) )  %>%
        as.data.frame() #%>%
      #add_column(population = pops, .after = 1)
    }
    #
    if( clim_var == "airt"){
      clim_m   <- clim_d %>%
        group_by(year, month) %>%
        summarise( ppt = mean(ppt, na.rm=T) )  %>%
        as.data.frame()
    }  
    # throw error
    if( !any( clim_var %in% c("precip","pet","airt")) ) {
      stop( paste0(clim_var," is not a supported varible") ) 
    }
    
    return(clim_m)
    
  }
  
  # monthly climate
  clim_m_l    <- lapply(clim_separate, clim_month)
  clim_m_m_l  <- lapply(clim_m_l, function(x) x$ppt)
  clim_m      <- do.call(cbind, clim_m_m_l)
  
  # get correlation of climates across sites 
  cor_mat           <- cor(clim_m) 
  
  cor_mat[ upper.tri(cor_mat) ]
  
}

# spp with multiple site information
mpl_spp_l   <- lapply(1:length(spp_list), spp_mpl_sites )
mpl_spp_df  <- Reduce(function(...) rbind(...), mpl_spp_l) %>%
                  subset( site_rep > 1 )
# site spp_diff
mpl_site_spp<- mpl_spp_df$species %>% as.character


# site level climates for all species 
site_clim_l     <- lapply(1:length(spp_list), spp_by_site_climate) 
site_clim_df    <- Reduce(function(...) bind_rows(...), site_clim_l) 

# between species differences in climate
spp_clims       <- site_clim_df %>% 
                      group_by( species ) %>%
                      summarise( mean_ppt = mean( mean_p ) )
spp_clim_diff   <- diff(spp_clims$mean_ppt) %>% as.matrix


# between site differences in climate
bw_site_diff_l  <- lapply(1:length(mpl_site_spp), bw_site_clim)
bw_site_diff    <- Reduce(function(...) rbind(...), bw_site_diff_l) 

# exclude species whose difference in climate is 0 (not a problem)
probl_spp <- bw_site_diff %>% 
                group_by( species) %>% 
                summarise( mean_spp_clim = mean(site_clim_diff) ) %>%
                ungroup %>%
                subset( mean_spp_clim != 0 ) %>%
                as.data.frame %>%
                .[,"species"] %>%
                as.character

# Cypripedium_fasciculatum and Brassica_insularis definitely raise a brow!
# Astragalus_scaphoides_2, Trillium_ovatum, 
subset(bw_site_diff, species %in% probl_spp )


# correlation of yearly climates
clim_corr_l <- lapply(mpl_site_spp, site_yr_clim_corr) %>%
                    setNames(mpl_site_spp)

# correlation of monthly climates
clim_m_corr_l <- lapply(mpl_site_spp, site_mon_clim_corr) %>%
                    setNames(mpl_site_spp)



# eucledian geographic distance --------------------------------------------------------------
euclid_dist <- function(ii){
  
  spp_name    <- probl_spp[ii]
  
  # site location
  sites_loc   <- lambdas %>% 
                    subset( SpeciesAuthor == spp_name ) %>%
                    dplyr::select(SpeciesAuthor,MatrixPopulation, Lat, Lon) %>%
                    unique 
    
  # site distance in kms
  site_km_dist<- sites_loc %>%
                    dplyr::select(Lat,Lon) %>%
                    as.matrix %>%
                    sp::spDistsN1( .[1,], longlat=TRUE )

  return( data.frame( species      = spp_name,
                      site_dist_km = site_km_dist)
        )

}
within_site_diff_l  <- lapply(1:length(probl_spp), euclid_dist)
within_site_diff    <- Reduce( function(...) rbind(...), within_site_diff_l)

# visualize species with site differences above 20/30Km
# species with site differences above 20/30 Km
spp_check <- within_site_diff %>% subset( site_dist_km > 20) %>% .$species %>% unique

# species Lat/Lon
spp_lat_lon <- lambdas %>%
                  subset( SpeciesAuthor %in% spp_check ) %>%
                  dplyr::select( SpeciesAuthor, Lat, Lon) %>%
                  unique


# vision species one by one
library(leaflet)
library(mapview)

# Astragalus_scaphoides_2: might show different climate - far away across mountains
leaflet(data = subset(spp_lat_lon, SpeciesAuthor == spp_check[1]) ) %>% 
  addTiles() %>% 
  leaflet::addCircleMarkers(~Lon, ~Lat)

# # 'leaflet' objects (image above)
# setwd("C:/cloud/MEGA/Projects/sApropos")
# mapshot(m, file = paste0(spp_check[1],".png"))

# Astragalus_tyghensis (all in eastern Oreon plains, but close to mountains)
leaflet(data = subset(spp_lat_lon, SpeciesAuthor == spp_check[2]) ) %>% 
  addTiles() %>% 
  leaflet::addCircleMarkers(~Lon, ~Lat)

# Brassica_insularis: Too different positions - All across Corsica!!! 
leaflet(data = subset(spp_lat_lon, SpeciesAuthor == spp_check[3]) ) %>% 
  addTiles() %>% 
  leaflet::addCircleMarkers(~Lon, ~Lat)

# Cirsium_pitcheri_4: Across great lake. No big differences IMO 
leaflet(data = subset(spp_lat_lon, SpeciesAuthor == spp_check[4]) ) %>% 
  addTiles() %>% 
  leaflet::addCircleMarkers(~Lon, ~Lat)

# Cypripedium_fasciculatum: Too different, from cascades to almost coastal mountains 
leaflet(data = subset(spp_lat_lon, SpeciesAuthor == spp_check[5]) ) %>% 
  addTiles() %>% 
  leaflet::addCircleMarkers(~Lon, ~Lat)

# Trillium_ovatum: Should experience diff. climate, around Missoula, so large topographic differences.
leaflet(data = subset(spp_lat_lon, SpeciesAuthor == spp_check[6]) ) %>% 
  addTiles() %>% 
  leaflet::addCircleMarkers(~Lon, ~Lat)
