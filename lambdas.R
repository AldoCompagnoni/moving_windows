# subsetting studies in the COMPADRE with enough temporal replication 
setwd("~/Dropbox/sAPROPOS project/DemogData")
library(dplyr)
library(tidyr)
library(testthat)
load("COMPADRE_v.X.X.X.2.Rdata")


# Select species with at least "tn" transitions ---------------------------------------

# number of transition 
tn <- 6

# COMPADRE ids of unmanipulated species with study duration > tn
spp_id  <- (compadre$metadata$StudyDuration >= tn &         # tn years of data
            compadre$metadata$MatrixTreatment == "Unmanipulated"  &
            compadre$metadata$MatrixComposite == "Individual" 
                # logical conditions to work on vital rate data
            # compadre$metadata$SurvivalIssue <= 1 &
            
            # compadre$metadata$MatrixFec == "Yes" &
            # compadre$metadata$MatrixSplit == "Divided" &
            ) %>% which 


# get lambda values ------------------------------------------------------------

# extract asymptotic lambda
lam_by_id <- function(ids, db_mat){
  
  matA <- db_mat[ids][[1]]$matA
  if( any( is.na(matA) ) ){
    return(NA)
  } else{
    return( Re(eigen(matA)$value[1]) )  
  }
  
}

# set up data frame: metadata + lambdas 
lambdas <- data.frame( lambda = sapply(spp_id, lam_by_id, compadre$mat) )
lam_df  <- bind_cols(compadre$metadata[spp_id,], lambdas) %>%
            subset(!is.na(lambda) ) %>% # remove NAs in lambda          
            arrange(SpeciesAuthor, MatrixPopulation, MatrixEndYear) 


# keep only species with at least **tn** lambdas per EACH MaxtrixPopulation -------------------------------
tn_yrs <- function(x){ sum( diff(x) == 1) >= (tn-1) }

# species and populations consecutive yearly data
consec_rep <- lam_df %>%
                group_by(SpeciesAuthor, MatrixPopulation) %>%
                summarise( consec_yr = tn_yrs(MatrixEndYear) ) %>%
                subset(consec_yr == TRUE) %>%
                select(SpeciesAuthor, MatrixPopulation) %>%
                unique

# select species X population with at least **tn** lambdas
lam_kp  <- subset(lam_df, SpeciesAuthor %in% consec_rep$SpeciesAuthor & 
                          MatrixPopulation %in% consec_rep$MatrixPopulation)

# tests ----------------------------------------------------------------------

# replication by species
rep_tst <- lam_kp %>%
            group_by(SpeciesAuthor, MatrixPopulation) %>%
            summarise( n = n() )
expect_equal(min(rep_tst$n), tn)

# are the *tn* years consecutive? 
consec_yr <- function(r, lam_kp, consec_rep){
  
  subset(lam_kp, SpeciesAuthor %in% consec_rep$SpeciesAuthor[r] & 
               MatrixPopulation %in% consec_rep$MatrixPopulation[r]) %>%
              .$MatrixEndYear %>%
              tn_yrs
  
}
cons_tst <- sapply(1:nrow(consec_rep), consec_yr, lam_kp, consec_rep)
expect_true( all(cons_tst) )


# store -----------------------------------------------------------------------------------
write.csv(lam_kp, paste0("lambdas_",tn,"tr.csv"), row.names = F)
