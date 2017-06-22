setwd("C:/Users/Aldo/Dropbox/sAPROPOS project/")
library(climwin)
library(dplyr)
library(tidyr)
library(rjags)

# read data
d       <- read.csv("DemogData/lambdas.csv")
precip  <- read.csv("DemogData/climate_data/monthly_ppt_Dalgleish.csv")

# analysis parameters  
spp     <- 4
m_back  <- 24

# format species ---------------------------------------------------

# select species
spp_dur <- d %>% 
            group_by(SpeciesAccepted) %>% 
            summarise(duration = length(unique(MatrixStartYear)))

# species list
spp_list<- d %>%
            subset(Lat == 38.8 & Lon == -99.2) %>%
            .[,"SpeciesAccepted"] %>%
            as.character %>%
            unique

# lambdas
xx     <- d %>%
            subset(SpeciesAccepted == spp_list[spp]) %>%
            dplyr::select(MatrixEndYear, lambda) %>%
            setNames(c("year","lambda")) %>%
            mutate(log_lambda = log(lambda))


# format climate ---------------------------------------------------------

# detrend climate
m_means <- colMeans(precip, na.rm=T)[-1]
d_prec  <- apply(precip[,-1], 2, FUN = scale, center = T, scale = T) 
det_prec<- cbind(precip[,1],d_prec) %>% as.data.frame
names(det_prec)[1] <- "YEAR"           

# select precipitation range
pr_long <- det_prec %>%
            subset(YEAR < 1974 & YEAR > 1934) %>%
            gather(month,precip,JAN:DEC) %>%
            setNames(c("year", "month", "precip")) %>%
            mutate(month = factor(month, levels = toupper(month.abb)) ) %>%
            mutate(month = as.numeric(month)) %>%
            arrange(year, month)

# array for # months before each sampling date.
prec_form <- function(x,dat,var){
  
  id <- which(dat$year == x & dat$month == 7) #"JUL"
  r  <- c(id:(id-(m_back-1)))
  return(dat[r,var])
  
}

# calculate monthly precipitation values
years   <- xx$year %>% unique
years   <- years[order(years)]
prec_l  <- lapply(years, prec_form, pr_long, "precip")

# put all in matrix form 
mat_form<- function(x) { 
  do.call(cbind, x) %>% 
    as.data.frame %>%
    setNames(years) 
}
x_prec  <- mat_form(prec_l)
x_prec  <- x_prec / diff(range(x_prec))

# Parametric climate weights --------------------------------------------
sink("Code/climwin/COMPADRE/norm_COMPADRE.txt")
cat("
    model{
      
      # priors for weights of monthly climate 
      m_c   ~ dunif(0.1, 24) 
      s_c   ~ dunif(0.1, 10) # Residual sd growth    
      tau_c <- 1 / (s_c * s_c) 
      
      # regression priors
      alpha ~ dunif(-10, 10) # Prior: intercept
      beta_c ~ dunif(-10, 10) # Prior: slope
      sigma ~ dunif(0, 10) # Residual sd     
      tau <- 1 / (sigma * sigma) 
      
      # Likelihoods
      
      # month weights 
      for(i in 1:W){
        #ws_raw[i] <- dweib(int[i],m_c,s_c)
        ws_raw[i] <- dnorm(int[i],m_c,tau_c)
      }

      # weights sum to 1
      ws_sum <- sum(ws_raw)
      ws     <- ws_raw / ws_sum

      # effect of climate on yearly population growth rates
      for(yr in 1:Y){ # loop over years
        # predictor as weighted climate predictors
        x_s[yr]  <- sum(x_prec[,yr] * ws)
      }      
      
      # loop through multiple data within each year
      for(n in 1:N){
        mu[n] <- alpha + beta_c*x_s[yr_i[n]]
        y1[n] ~ dnorm(mu[n],tau)
      }
      
    }#End model",fill=T)
sink()


# bundle data
# index for year
yr_i      <- xx$year %>% as.factor %>% as.numeric
jag_data  <-list(y1 = xx$log_lambda,
                 x_prec = x_prec[1:m_back,],
                 int = c(1:m_back),
                 N = nrow(xx),
                 Y = ncol(x_prec),
                 yr_i = yr_i,
                 W = nrow(x_prec))

# Parameters to estimate
parameters <- c("m_c", "s_c",
                "alpha",  "beta_c","sigma")

#Initial values
inits <- function(){list( m_c = runif(1, 0.11, m_back),#runif(1, 0.11, m_back), 
                          s_c = rlnorm(1, 0, sdlog = 0.01),
                          alpha = rnorm(1, 0, 2),
                          beta_c  = rnorm(1, 0, 2), 
                          sigma = rlnorm(1)) }

## MCMC settings
ni<-100000
nb<-50000
nt<-5
nc<-3

# JAG model
mod <- jags.model(file="Code/climwin/COMPADRE/norm_COMPADRE.txt", data=jag_data, 
                  inits=inits, n.chains = nc, quiet=F)

# "fit" model
out       <- coda.samples(mod, variable.names=parameters, 
                          n.iter=ni, thin = nt)
# keep samples only after burnin
out_keep  <- window(out, start=nb, end=ni)
plot(out_keep)


# plot ---------------------------------------------------------------

# extract posterior samples
sim_extract <- function(i) {
  out_keep[[i]] %>% 
    as.data.frame
}

sim_l     <- lapply(c(1:nc), sim_extract)
sim_df    <- do.call(rbind,sim_l)
pars_id   <- seq(1, nrow(sim_df), length.out = 100) %>%
              round()
pars_weib <- sim_df[pars_id,]
sim_means <- apply(sim_df,2,mean)


# figure: sensitivity biomass ~ monthly precip --------------------------------
tiff(paste0("Code/climwin/COMPADRE/results/",gsub(" ","_",spp_list[spp]),"_vs_precip.tiff"),
     unit="in", width=6.3,height=3.5,res=600,compression="lzw")

# sequence of months (from 1st all the way to "back")
month_x <- seq(1,m_back,by=0.1)

# store 100 normal distrib. across 1:m_back months
store_norm <- function(i){
  dnorm(month_x, pars_weib[i,"m_c"],pars_weib[i,"s_c"])
}
# store 100 posterior samples of month sensitivity
month_sens <- do.call(rbind, lapply(c(1:100), store_norm) )

# month sensitivities
mean_m_sen<- dnorm(seq(1,m_back,1), sim_means["m_c"], sim_means["s_c"])
# climate predictor
clim_mat   <- sweep(x_prec, 1, mean_m_sen, "*") %>%
                apply(2,sum)
# plotting
clim_x      <- stack(clim_mat) %>%
                setNames(c("clim","year")) %>%
                mutate(year = as.numeric(as.character(year)))
plot_mat    <- merge(xx, clim_x)

# plot 100 lines
abline(a = sim_means["alpha"],
       b = sim_means["beta_c"], lwd = 2)



# plots ------------------------------------------------------------------
par(mfrow=c(1,2), mar=c(3.2,3.2,0.3,0.1), mgp = c(2,0.7,0), cex.lab = 1.2)
plot(month_x,
     seq(1,m_back,length.out=length(seq(1,m_back,by=0.1))),
     type="n", ylim = c(0,max(month_sens)),
     xlab = "Month", ylab = "Sensitivity to monthly precip.")

# plot 100 posterior samples
lapply(c(1:100), function(i) 
                 lines(month_x, 
                       month_sens[i,], 
                       lwd=0.5, col = "grey") )

# plot mean
lines(month_x,
      dnorm(month_x, sim_means["m_c"], sim_means["s_c"]),
      type="l", lwd=2)

# lambda ~ climate predictor 
plot(log_lambda ~ clim, data = plot_mat, pch = 16,
     xlab = "Precip. (weighted avg.)",
     ylab = expression(lambda))
# posterio of predictors 
lapply(1:100, function(i) 
              abline(a = pars_weib[i,"alpha"], col = "grey",
                     b = pars_weib[i,"beta_c"], lwd = 2)
       )
# best model fit
abline(a = sim_means["alpha"],
       b = sim_means["beta_c"], lwd = 2)
points(log_lambda ~ clim, data = plot_mat, pch = 16,
     xlab = "Precip. (weighted avg.)",
     ylab = expression(lambda), new=T)

dev.off()

