# subsetting studies in the COMPADRE with enough temporal replication 
setwd("~/Dropbox/sAPROPOS project/DemogData")
library(dplyr)
library(testthat)
load("COMPADRE_v.X.X.X.2.Rdata")


# Select species with at least "tn" transitions ---------------------------------------

# number of transition 
tn <- 6

# COMPADRE ids of unmanipulated species with study duration > tn
spp_id  <- (compadre$metadata$StudyDuration >= tn &         # tn years of data
            compadre$metadata$MatrixTreatment == "Unmanipulated"
                # logical conditions to work on vital rate data
            # compadre$metadata$SurvivalIssue <= 1 &
            # compadre$metadata$MatrixComposite == "Individual"  &
            # compadre$metadata$MatrixFec == "Yes" &
            # compadre$metadata$MatrixSplit == "Divided" &
            ) %>% which 


# get lambda values ------------------------------------------------------------

# extract asymptotic lambda
lam_by_id <- function(ids, db_mat){
  
  matA <- db_mat[ids][[1]]$matA
  if(any( is.na(matA) ) ){
    return(NA)
  } else{
    return( Re(eigen(matA)$value[1]) )  
  }
  
}

# set up data frame: metadata + lambdas 
lambdas <- data.frame( lambda = sapply(spp_id, lam_by_id, compadre$mat) )
lam_df  <- bind_cols(compadre$metadata[spp_id,], lambdas) %>%
            subset(!is.na(lambda) ) # remove NAs in lambda


# keep only species with at least **tn** lambdas -------------------------------

# names of species to keep
spp_nam <- lam_df %>%
            group_by(SpeciesAuthor) %>%
            summarise( n = n() ) %>%
            subset(n >= tn) %>%
            select(SpeciesAuthor) %>%
            .$SpeciesAuthor %>% unique

# select only species with at least **tn** lambdas
lam_kp  <- subset(lam_df, SpeciesAuthor %in% spp_nam)

# test replication by species
rep_tst <- lam_kp %>%
            group_by(SpeciesAuthor) %>%
            summarise( n = n() )
expect_equal(min(rep_tst$n), tn)


# out 
write.csv(lam_kp, paste0("lambdas_",tn,"tr.csv"), row.names = F)
