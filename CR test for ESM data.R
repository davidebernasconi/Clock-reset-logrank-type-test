###
#### my_ESM_boot_logrank: "Clock-reset logrank-type" test for extended semi-markov data ####
###

library(boot)
my_ESM_boot_logrank<-function(km.fit,clockbacktime,eventindicator,expbetax) { #arguments: km.fit, clockbacktime, eventindicator e expbetax 
  
  U.func <- function(x,i) {
    htdata<-x[i,]
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
    m<-coxph(Surv(Tdelta,Death.tv>0)~T.htx,data=htdatasplitHTx[htdatasplitHTx$HTx.tv==1,],x=T)
    clockbacktime<-htdatasplitHTx$Tdelta[htdatasplitHTx$HTx.tv==1]
    eventindicator<-htdatasplitHTx$Death.tv[htdatasplitHTx$HTx.tv==1]
    expbetax<-exp(m$coefficients*htdatasplitHTx$T.htx[htdatasplitHTx$HTx.tv==1])
    x<-km.fit
    event.timesA<-ifelse(x$n.event==0,NA,x$time)
    event.timesB<-unique(clockbacktime[eventindicator>0])
    event.times<-c(event.timesA,event.timesB)
    event.times<-sort(event.times)
    event.times<-unique(event.times)
    vect.U<-NULL
    for (j in event.times) {
      if (j>max(x$time) | j>max(clockbacktime[eventindicator>0])) {} else {
        nA<-ifelse(j<min(x$time),x$n.risk[x$time==min(x$time)],ifelse(j %in% x$time,x$n.risk[x$time==j],x$n.risk[x$time==min(x$time[x$time>j])]))
        nB<-sum(expbetax[clockbacktime>=j])
        dA<-ifelse(j %in% x$time,x$n.event[x$time==j],0)
        dB<-sum(eventindicator[clockbacktime==j])
        d<-dA+dB
        n<-nA+nB 
        Expected.dA<-nA*(d/n)
        U.temp<-dA-Expected.dA
        vect.U<-c(vect.U,U.temp)
      }
    }
    U<-sum(vect.U)
    return(U)  
  }
  U.boot <- boot(htdata, U.func, R=1000)
  
  var.U<-var(U.boot$t)
  U<-U.boot$t0
  Q<-(U^2)/var.U
  p.value<-pchisq(Q, df=1, ncp = 0, lower.tail = F, log.p = FALSE)
  output<-list(Q=Q,p.value=p.value)
  return(output)
}
