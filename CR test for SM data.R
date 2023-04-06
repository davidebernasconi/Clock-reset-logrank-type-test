###
#### my_SM_logrank: "Clock-reset logrank-type" test for semi-markov data ####
###

my_SM_logrank<-function(km.fit,clockbacktime,eventindicator,expbetax) { #arguments: km.fit, clockbacktime, eventindicator e expbetax 
  x<-km.fit
  event.timesA<-ifelse(x$n.event==0,NA,x$time)
  event.timesB<-unique(clockbacktime[eventindicator>0])
  event.times<-c(event.timesA,event.timesB)
  event.times<-sort(event.times)
  event.times<-unique(event.times)
  vect.U<-NULL
  vect.var.U<-NULL
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
      var.dA.temp<-(nA*(d/n)*(1-d/n))*((n-nA)/(n-1))
      vect.U<-c(vect.U,U.temp)
      vect.var.U<-c(vect.var.U,var.dA.temp)
    }
  }
  U<-sum(vect.U)
  var.U<-sum(vect.var.U)
  Q<-(U^2)/var.U
  p.value<-pchisq(Q, df=1, ncp = 0, lower.tail = F, log.p = FALSE)
  output<-list(Q=Q,p.value=p.value)
  return(output)
}