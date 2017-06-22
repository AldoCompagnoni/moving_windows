bjtwd<-"C:/Users/admin_bjt162/Dropbox/A.Current/Ongoing_Collab_Research/sApropos project/"
#afcwd<-"C:/Users/Aldo/Dropbox/sAPROPOS project/"
setwd(bjtwd)
spp_list

library(climwin)
library(dplyr)
library(tidyr)
library(rjags)
library(dismo)

# read data
d       <- read.csv("DemogData/lambdas.csv")
precip  <- read.csv("DemogData/climate_data/monthly_ppt_Dalgleish.csv")

for(ii in c(1,3)){

  # analysis parameters  
  spp     <- ii
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
  
  pmat<-t(x_prec)
  
  lags <- matrix(0,nrow(pmat),ncol(pmat)); 
  for(i in 1:ncol(lags)) lags[,i]=i; 
  lags=as.matrix(lags); 
  
  pmat <- pmat/diff(range(pmat));
  
  dat<-cbind(xx,pmat)
  dat<-as.data.frame(dat, colnames=TRUE)
  dat$lags<-lags
  dat$pmat<-pmat
  
  mod2<-gam(log_lambda ~ s(lags, k=10, by=pmat, bs="cs"),
            data=dat,method="GCV.Cp",gamma=1.4) 
  
  tiff(file=paste0("Code/climwin/COMPADRE/results/",spp_list[ii],"_spline.tiff"), 
       unit="in", width=12,height=7,res=400,compression="lzw")
  
  par(mfrow=c(1,2))
  plot(mod2,main=spp_list[ii])
  abline(h=0)
  
  uniYears<-rownames(pmat)
  predNull<-rep(NA,times=length(uniYears))
  pred2<-rep(NA,times=length(uniYears))
  
  for (i in 1:nrow(pmat)){
    dati<-dat[-i,]
    predNull[i]<-mean(dat$log_lambda[-i])
    
    mod2.loo<-gam(log_lambda ~ s(lags, k=10, by=pmat, bs="cs",sp=mod2$sp),
                  method="GCV.Cp",gamma=1.4, data=dati) 
    
    pred2[i]<-predict(mod2.loo,newdata = dat[i,], type="response")
    print(c(ii,i))
  }

  winPred<-read.csv(paste0("Code/climwin/COMPADRE/results/moving_window_spp",ii,"_pred.csv"),header=T)
  
  dev0<-calc.deviance(predNull, dat$log_lambda, 
                      weights = rep(1,length(pred2)),  
                      family="gaussian", calc.mean = TRUE)  
  dev1<-calc.deviance(pred2, dat$log_lambda, 
                weights = rep(1,length(pred2)),  
                family="gaussian", calc.mean = TRUE)
  devg<-calc.deviance(winPred$pred_gaus, dat$log_lambda, 
                      weights = rep(1,length(pred2)),  
                      family="gaussian", calc.mean = TRUE)
  deve<-calc.deviance(winPred$pred_expp, dat$log_lambda, 
                      weights = rep(1,length(pred2)),  
                      family="gaussian", calc.mean = TRUE)
  devi<-calc.deviance(winPred$pred_c1, dat$log_lambda, 
                      weights = rep(1,length(pred2)),  
                      family="gaussian", calc.mean = TRUE)
  dev24<-calc.deviance(winPred$pred_c2, dat$log_lambda, 
                      weights = rep(1,length(pred2)),  
                      family="gaussian", calc.mean = TRUE)
  
  plot(dat$log_lambda~as.numeric(as.character(uniYears)), type="o", ylab="logLambda", xlab="Left-out year", col="black", main=spp_list[ii], ylim=c(-3,4))
  points(predNull~as.numeric(as.character(uniYears)), type="o", col="red")
  points(winPred$pred_c1~as.numeric(as.character(uniYears)), type="o", col="yellow")
  points(winPred$pred_c2~as.numeric(as.character(uniYears)), type="o", col="orange")
  points(pred2~as.numeric(as.character(uniYears)), type="o", col="blue")
  points(winPred$pred_expp~as.numeric(as.character(uniYears)), type="o", col="green")
  points(winPred$pred_gaus~as.numeric(as.character(uniYears)), type="o", col="purple")

  dev.0=round(dev0,3)
  dev.1=round(dev1,3)
  dev.g=round(devg,3)
  dev.e=round(deve,3)
  dev.i=round(devi,3)
  dev.24=round(dev24,3)
  
  legend("topright", pch=1,lty=1, col=c("black","blue","red","purple","green"),
         legend=c("Real loglambda", paste("Climate Spline; dev=",dev.1),paste("Nothing; dev=",dev.0), paste("Window G; dev=",dev.g),
                  paste("Window E; dev=",dev.e),paste("Window I; dev=",dev.i),paste("Window 24; dev=",dev.24)), bty="n")
  
  dev.off()
}
