#### set wd ####
#setwd("C:\\...")

#### call the functions ####
source("CR test for SM data.R")
source("CR test for ESM data.R")

#### read data ####

htdata<-read.table("htdata.txt",header=T,sep="")

# select only subjects with status=2B (lowest priority to transplant)
htdata<-htdata[htdata$status=="2B",]

table(htdata$HTx)
summary(htdata$T.death[htdata$HTx==1])
summary(htdata$T.htx[htdata$HTx==1])

table(htdata$HTx,htdata$Death)
summary(htdata$T.death)

#check Markov
library(survival)
m<-coxph(Surv(T.htx,T.death,Death)~T.htx,data=htdata[htdata$HTx==1,])
summary(m)
summary(m)$sctest

#check semi-Markov
library(survival)
m<-coxph(Surv(T.death-T.htx,Death)~T.htx,data=htdata[htdata$HTx==1,])
summary(m)
summary(m)$sctest

# check assumptions
library(survminer)
ggcoxzph(cox.zph(m),ylab="Beta(t)")

library(Greg)
library(splines)

m<-coxph(Surv(T.death-T.htx,Death)~bs(T.htx,3),data=htdata[htdata$HTx==1,])
par(mar=c(4,4,1,1),mfrow=c(1,1))
plotHR(m, term = "T.htx", plot.bty = "o", xlim = c(0, 3), xlab = "Waiting time to HTx")


#### rearrange data and apply CR and MB tests ####

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


# MB-test
m<-coxph(Surv(Tstart,Tstop,Death.tv>0)~htdatasplitHTx$HTx.tv,data=htdatasplitHTx,timefix = FALSE)
summary(m)
summary(m)$sctest



# CR-test
km.fit<-survfit(Surv(Tstop,Death.tv>0)~1,data=htdatasplitHTx[htdatasplitHTx$HTx.tv==0,])
m<-coxph(Surv(Tdelta,Death.tv>0)~T.htx,data=htdatasplitHTx[htdatasplitHTx$HTx.tv==1,],x=T)
clockbacktime<-htdatasplitHTx$Tdelta[htdatasplitHTx$HTx.tv==1]
eventindicator<-htdatasplitHTx$Death.tv[htdatasplitHTx$HTx.tv==1]
#expbetax<-exp(m$coefficients*htdatasplitHTx$T.htx[htdatasplitHTx$HTx.tv==1])
expbetax<-exp(0*htdatasplitHTx$T.htx[htdatasplitHTx$HTx.tv==1])

LR<-my_SM_logrank(km.fit,clockbacktime,eventindicator,expbetax)
LR

#ECR-test
km.fit<-survfit(Surv(Tstop,Death.tv>0)~1,data=htdatasplitHTx[htdatasplitHTx$HTx.tv==0,])
m<-coxph(Surv(Tdelta,Death.tv>0)~T.htx,data=htdatasplitHTx[htdatasplitHTx$HTx.tv==1,],x=T)
clockbacktime<-htdatasplitHTx$Tdelta[htdatasplitHTx$HTx.tv==1]
eventindicator<-htdatasplitHTx$Death.tv[htdatasplitHTx$HTx.tv==1]
expbetax<-exp(m$coefficients*htdatasplitHTx$T.htx[htdatasplitHTx$HTx.tv==1])

set.seed(7)
LR<-my_ESM_boot_logrank(km.fit,clockbacktime,eventindicator,expbetax)
LR

#### plot hazards by quartiles of waiting time vs hazard without HTx ####
library(bshazard)
summary(htdatasplitHTx$T.htx[htdatasplitHTx$HTx.tv==1])
htdatasplitHTx$groups<-ifelse(htdatasplitHTx$HTx.tv==0,NA,
                              ifelse(htdatasplitHTx$T.htx<=quantile(htdatasplitHTx$T.htx[htdatasplitHTx$HTx.tv==1],0.25),1,
                                     ifelse(htdatasplitHTx$T.htx<=quantile(htdatasplitHTx$T.htx[htdatasplitHTx$HTx.tv==1],0.5),2,
                                            ifelse(htdatasplitHTx$T.htx<=quantile(htdatasplitHTx$T.htx[htdatasplitHTx$HTx.tv==1],0.75),3,4))))
table(htdatasplitHTx$groups)

fit.noHTx<-bshazard(Surv(Tstop,Death.tv) ~ 1,data=htdatasplitHTx[htdatasplitHTx$HTx.tv==0,])
par(mar=c(4,4,1,1))
plot(fit.noHTx,ylim=c(0,0.5),col="darkgrey",lwd=3,conf.int = F,ylab="Mortality rate (hazard)",xlab="Time (years)")

fit1<-bshazard(Surv(Tdelta,Death.tv) ~ 1,data=htdatasplitHTx[htdatasplitHTx$groups==1,])
lines(fit1,col=1,conf.int = F)
fit2<-bshazard(Surv(Tdelta,Death.tv) ~ 1,data=htdatasplitHTx[htdatasplitHTx$groups==2,])
lines(fit2,col=2,conf.int = F)
fit3<-bshazard(Surv(Tdelta,Death.tv) ~ 1,data=htdatasplitHTx[htdatasplitHTx$groups==3,])
lines(fit3,col=3,conf.int = F)
fit4<-bshazard(Surv(Tdelta,Death.tv) ~ 1,data=htdatasplitHTx[htdatasplitHTx$groups==4,])
lines(fit4,col=4,conf.int = F)
legend("topright",c("No HTx","HTx after (0-1.8] months","HTx after (1.8-5.4] months",
                    "HTx after (5.4-12] months","HTx after >12 months"),
       col=c("darkgrey",1,2,3,4),lwd=c(3,1,1,1,1))
