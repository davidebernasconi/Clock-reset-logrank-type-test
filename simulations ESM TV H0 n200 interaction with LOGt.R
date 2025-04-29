### 
### EXTENDED SEMI MARKOV DATA (time-varying effect of waiting time)
### 

### !!! ECR test with TV effect of waiting time (LOGt) !!! ###
 
source("LOGt CR test for ESM TV data.R")

#----------------------------------------------#
####  simulations under the null H0, n=200  ####
#----------------------------------------------#

# hazard of death (and potential surv) before HTx and after Htx at time 0 are the same

### parameters set and true quantities
t<- seq(0, 10, 0.01)
# exponential hazard
g<-0.5
haz.list.htx<-g*t^0
surv.list.htx<-exp(-g*t)
# weibull hazards
k<-0.6
q<-0.7
b<-(-0.2)
haz.list.death<-q*k*((k*t)^(q-1))
surv.list.death<-exp(-((k*t)^q))


### data generation

library(survival)

set.seed(123456789)
nsim<-1000
p.value<-NULL
for (i in 1:nsim) {
  # set number of observations
  n<-200
  
  # generate latent censoring times
  t.cens<- runif(n, min=0.25, max=8)
  
  # generate latent times from list entry to transplant
  U<-runif(n, min=0, max=1)
  t.htx<- ((-log(U))/g)
  # generate latent times from list entry to death
  U<-runif(n, min=0, max=1)
  t.list.death<- (1/(k))*((-log(U)))^(1/q)
  # generate latent times from transplant to death
  U<-runif(n, min=0, max=1)
  t.htx.death<- ( (((-log(U))*(b*t.htx+1))^(1/q))/((5^(-b*t.htx))*k) )^(1/(b*t.htx+1))  
  t.htx.death<-ifelse(is.na(t.htx.death),0,t.htx.death) 
  # correzione necessaria perchè alcuni tempi molti piccoli (vicini a 0) vengono settati a NA...
  
  # observed time and event indicators
  htx<-ifelse(t.htx<t.list.death, 1, 0)
  t.death<-ifelse(htx==0, t.list.death, t.htx+t.htx.death)
  T.death<-pmin(t.death,t.cens)
  Death<-ifelse(t.death<t.cens,1,0)
  HTx<-ifelse(T.death<t.htx,0,htx)
  T.htx<-ifelse(HTx==0,NA,pmin(t.htx,T.death))
  T.htx.death<-ifelse(HTx==0,NA,ifelse(HTx==1 & Death==0,T.death-t.htx,t.htx.death))
  
  # put all in a dataset
  htdata<-as.data.frame(cbind(T.death,Death,T.htx,T.htx.death,HTx,id=seq(1,n)))
  
  
  ### data analysis
  
  # rearrange data: split follow-up time of transplanted in two parts: pre and post HTx
  htdata_HTx<-htdata[htdata$HTx==1,]
  htdata_HTx<-htdata_HTx[rep(row.names(htdata_HTx), 2),]
  htdatasplitHTx<- rbind(htdata[htdata$HTx==0,],htdata_HTx)
  htdatasplitHTx$split<-c(rep("pre HTx",dim(htdata)[1]),rep("post HTx",dim(htdata[htdata$HTx==1,])[1]))
  # times in "counting process" format
  htdatasplitHTx$Tstart<-ifelse(htdatasplitHTx$split=="pre HTx",0,htdatasplitHTx$T.htx)
  htdatasplitHTx$Tstop<-ifelse(htdatasplitHTx$HTx==1 & htdatasplitHTx$split=="pre HTx",htdatasplitHTx$T.htx,htdatasplitHTx$T.death)
  htdatasplitHTx$Tdelta<-htdatasplitHTx$Tstop-htdatasplitHTx$Tstart
  # time-varying indicators for transplanted
  htdatasplitHTx$HTx.tv<-ifelse(htdatasplitHTx$split=="pre HTx" & htdatasplitHTx$HTx==1,0,htdatasplitHTx$HTx)
  htdatasplitHTx$Death.tv<-ifelse(htdatasplitHTx$split=="pre HTx" & htdatasplitHTx$HTx==1,0,htdatasplitHTx$Death)
  
  km.fit<-survfit(Surv(Tstop,Death.tv>0)~1,data=htdatasplitHTx[htdatasplitHTx$HTx.tv==0,])
  
  htdatasplitHTx$Tdelta<-ifelse(htdatasplitHTx$Tdelta==0,htdatasplitHTx$Tdelta + 0.00000001,htdatasplitHTx$Tdelta)
  htdatasplitHTx.1<-htdatasplitHTx[htdatasplitHTx$HTx.tv==1,]
  htdatasplitHTx.1$idnew<-1:dim(htdatasplitHTx.1)[1]
  htdataHTX.TV <- survSplit(Surv(Tdelta,Death.tv) ~., htdatasplitHTx.1,
                            cut=unique(htdatasplitHTx.1$Tdelta))
  htdataHTX.TV$interact<-htdataHTX.TV$T.htx*log(htdataHTX.TV$Tdelta)
  m <- coxph(Surv(tstart, Tdelta, Death.tv>0) ~ T.htx + interact, data= htdataHTX.TV,x=T,
             control = coxph.control(timefix = FALSE))
  
  clockbacktime<-round(htdatasplitHTx$Tdelta[htdatasplitHTx$HTx.tv==1],8)
  eventindicator<-htdatasplitHTx$Death.tv[htdatasplitHTx$HTx.tv==1]
  
  expbetax<-exp(m$coefficients[1]*htdataHTX.TV$T.htx + m$coefficients[2]*htdataHTX.TV$interact )
  expbetax.df<-data.frame(idnew=htdataHTX.TV$idnew,time=htdataHTX.TV$Tdelta,expbetax)  
  
  LR<-my_ESM_boot_logrank_tv(km.fit,clockbacktime,eventindicator,expbetax.df)
  
  p.value<-c(p.value,LR$p.value)
  
  print(i)
}

summary(p.value)
table(p.value<0.05)
table(p.value<0.01)
table(p.value<0.1)
