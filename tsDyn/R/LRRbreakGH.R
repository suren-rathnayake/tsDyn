

##################
## Function
##################


LRRbreak.test<-function(y,x,model=c("C","CT","CS"), AdfLagChoice=c("User", "AIC", "BIC", "P-Val"),ADFlagmax=6, trim=0.15){ 
  #check args
  modelName<-match.arg(model)
  model<-switch(modelName, "C"=2, "CT"=3,"CS"=4)
  choice<-match.arg(AdfLagChoice)
  if(length(y)!=nrow(as.matrix(x)))    stop("X and y are not of the same size")
  
#set-up temporary variables
  n<-length(y)
  begin<-round(trim*n)
  final<-round((1-trim)*n)
  temp1<-c(rep(0,final-begin+1))
  temp4<-temp3<-temp2<-temp1

  t<-begin
  k<-ADFlagmax
##############
#ADF internal function
############
 
  
############
## Philips Za and Zt internal function
############
  GHphil <-function(yP,xP){
    nP<-nrow(as.matrix(xP))
    
#OLS regression 
    regP<-lm(yP~ -1 +xP)
    b<-coefficients(regP)
    eP<-residuals(regP)

# OLS regression on residuals
    regP2<-lm(eP[2:nP]~ -1 +eP[1:(nP-1)])			
    be<-coefficients(regP2)
    ue<-eP[2:nP]-eP[1:(nP-1)]*be

# calculate bandwidth number
    nu<-length(ue)	
    reg3<-lm(ue[2:nu]~ -1 + ue[1:(nu-1)])
    bu<-coefficients(reg3)
    uu<-ue[2:nu]-ue[1:(nu-1)]*bu
    su<-mean(uu^2)
    a2<-(4*(bu^2)*su/((1-bu)^8))/(su/(1-bu)^4)
    bandwidth<-1.3221*((a2*nu)^0.2)
    m<-bandwidth 
    j<-1
    lemda <-0
    
    while (j<=m){
      gama<-t(ue[1:(nu-j)])%*%ue[(j+1):nu]/nu
      c<-j/m
      w<-(75/(6*pi*c)^2)*(sin(1.2*pi*c)/(1.2*pi*c)-cos(1.2*pi*c))
      lemda<-lemda+w*gama
      j<-j+1
    }
    
    p<-sum(eP[1:(nP-1)]*eP[2:nP]-lemda)/sum(eP[1:(nP-1)]^2)
    za<-nP*(p-1)
    sigma2<-2*lemda+t(ue)%*%ue/nu
    s<-sigma2/(t(eP[1:(nP-1)])%*%eP[1:(nP-1)])
    zt<-(p-1)/sqrt(s)
    list(za=za, zt=zt)
  }
##############
##Compute ADF, Zt and Za for each t
##############
#   while (t<=final){
#   for(t in seq(begin, final)){
library(foreach)
temps<-foreach(t=seq(begin, final), .combine="rbind") if(hpc=="none") %do% else %dopar%{
    dummy<-c(rep(0, t), rep(1, n-t))
    x1<-switch(modelName, "C"=cbind(1, dummy,x), "CT"=cbind(1, dummy, 1:n, x), "CS"=cbind(1, dummy, x, x*dummy))
    #adjust regressors for different models
    if (model==3){ 
      x1<-cbind(rep(1,n), dummy, 1:n, x)}
    else if (model==4){
      x1<-cbind(rep(1,n), dummy, x, x*dummy)}
    else if (model==2){
      x1<-cbind(c(rep(1,n)), dummy,x)}
    
   #compute ADF for each t 
    adfTemp<-GHadf(y, x1, kmax=k,choice=choice)
    temp1<-adfTemp$tstat
    temp2<-adfTemp$lag
#     temp1[t-begin+1]<-adfTemp$tstat
#     temp2[t-begin+1]<-adfTemp$lag

#compute Phillips Zt and Za for each t 
    tempPP<-GHphil(y,x1)
    temp3<-tempPP$za
    temp4<-tempPP$zt
c(temp1, temp2, temp3, temp4)
#     temp3[t-begin+1]<-tempPP$za
#     temp4[t-begin+1]<-tempPP$zt

#     t<-t+1
}

# Min ADF test and breakpoint
  minADF<-min(temps[,1])
  lagminADF<-which.min(temps[,1])
  breakpointADF<-(lagminADF+begin-1)/n
  if(class(y)=="ts") breakpointADF<-time(y)[round(breakpointADF*n)]
  lag<-temps[lagminADF,2]	

# Min Phillips Za and Zt test and breakpoint
  minZa <-min(temps[,3])
  lagminZa<-which.min(temps[,3])
  breakpointZa <-(lagminZa+begin-1)/n
  if(class(y)=="ts")breakpointZa<-time(y)[round(breakpointZa*n)]

  minZt <-min(temps[,4])
  lagminZt<-which.min(temps[,4])
  breakpointZt<-(lagminZt+begin-1)/n
  if(class(y)=="ts")breakpointZt<-time(y)[round(breakpointZt*n)]

#Critical values for ADF and Zt, ordered following: m=1, C, C/T, C/S, m=2, C, C/T...
  table1<-rbind(
                c(5.13, 4.83, 4.61, 4.34, 2.25), c(5.45, 5.21, 4.99,4.72, 2.72), c(5.47, 5.19, 4.95,4.68,2.55),
c(5.44, 5.16, 4.92, 4.69, 2.61), c(5.80, 5.51, 5.29, 5.03, 3.01), c(5.97, 5.73, 5.50, 5.23, 3.12), 
c(5.77, 5.50, 5.28, 5.02, 2.96), c(6.05, 5.79, 5.57, 5.33, 3.33), c(6.51, 6.23, 6.00, 5.75, 3.65),
c(6.05, 5.80, 5.56, 5.31, 3.26), c(6.36, 6.07, 5.83, 5.59, 3.59), c(6.92,6.64, 6.41, 6.17, 4.12))

table1<- -table1
colnames(table1)<-c("0.01", "0.025", "0.05", "0.10", "0.975")

#Critical values for Za, ordered following: m=1, C, C/T, C/S, m=2, C, C/T...
  table2<-rbind(
c(50.07, 45.01, 40.48, 36.19, 10.63), c(57.28, 52.09, 47.96, 43.22, 15.90), c(57.17, 51.32, 47.04, 41.85, 13.15),
c(57.01, 51.41, 46.98, 42.49, 14.27), c(64.77, 58.57, 53.92, 48.94, 19.19), c(68.21, 63.28, 58.33, 52.85, 19.72), 
c(63.64, 57.96, 53.58, 48.65, 18.20), c(70.27, 64.26, 59.76, 54.94, 22.72), c(80.15, 73.91, 68.94, 63.42, 26.64),
c(70.18, 64.41, 59.40,54.38, 22.04), c(76.95, 70.56, 65.44, 60.12, 26.46),c(90.35, 84.00, 78.52, 72.56, 33.69))

  table2<- -table2
  colnames(table2)<-c("0.01", "0.025", "0.05", "0.10", "0.975")

  m<-ncol(as.matrix(x))
  if(m>4)m<-4

#Selection of the appropriate critical value for ADF and Zt
  selectCol<-switch(modelName, "C"=(m*3)-2, "CT"=(m*3)-1, "CS"=(m*3))
  crit.val1<-table1[selectCol,]
  crit.val2<-table2[selectCol,]
#   if(model==2)crit.val1<-table1[(m*3)-2,]
#   if(model==3)crit.val1<-table1[(m*3)-1,]
#   if(model==4)crit.val1<-table1[(m*3),]

#Selection of the appropriate critical value for Za
#   if(model==2)crit.val2<-table2[(m*3)-2,]
#   if(model==3)crit.val2<-table2[(m*3)-1,]
#   if(model==4)crit.val2<-table2[(m*3),]

 
#Output
  res<-list(adf=minADF, critical.values1=crit.val1, critical.values2=crit.val2, breakpointADF=breakpointADF, ADFlag=lag, 
Zt=minZt,breakpointZt=breakpointZt,Za=minZa,  breakpointZa=breakpointZa, model=modelName, vars=ncol(as.matrix(x)),
AdfLagChoice=AdfLagChoice, TestValues=cbind(temp1=temps[,1], temp3=temps[,3], temp4=temps[,4]))
  class(res) <- "GregHan1996"
  return(res)
}

##########
### ADF INTERNAL
### ########


 GHadf<-function(y,x, kmax, choice=c("User", "AIC", "BIC", "P-Val"), pvalue=0.05){
    n<-length(y)
    regadf1<-lm(y~ -1 +x)
    e<-residuals(regadf1)
    de<-e[2:n]-e[1:(n-1)] #difference of residuals
    ic<-0
    k<-kmax
    temp2<-temp1<-rep(0,kmax+1)
    
    while (k>=0){			
      yde<-de[(1+k):(n-1)]
      n1<-length(yde)
      xe<-cbind(e[(1+k):(n-1)], embed(de,k+1)[,-1])     # matrix for independent variable(lagged residuals)
      reg2<-lm(yde~ -1 +xe)
      tval <-summary(reg2)$coefficients[,3]
      pval <-summary(reg2)$coefficients[,4]	
      e1<-residuals(reg2)
      
	#lag selection
      if (choice=="User"){
        temp1[k+1]<- -1000
        temp2[k+1]<-tval[1]
        break}
      else if (choice=="AIC"){
        aic2<-AIC(reg2)
        aic<-log(t(e1)%*%e1)+2*((k+2)/n1)
        ic<-aic}
      else if (choice=="BIC"){
        bic<-log(t(e1)%*%e1)+(k+2)*log(n1)/n1
        bic2<-AIC(reg2, k=log(length(yde)))
        ic<-bic}
      else if (choice=="P-Val"){
        if (abs(pval[k+1])<=pvalue|k==0){
          temp1[k+1]<- -1000
          temp2[k+1]<-tval[1]
          break}
      }
      temp1[k+1]<-ic
      temp2[k+1]<-tval[1]
      k<-k-1	
    }
    lag<-which.min(temp1)
    tstat<-temp2[lag]
    lag1<-lag-1
    list(tstat=tstat, lag=lag1)
  }


#check GADF
if(FALSE){
  library(urca)
  GHadf(y,x, kmax=1)$tstat
  summary(ur.df(residuals(lm(y~x-1))))@teststat

  GHadf(y,x, kmax=2)$tstat
  summary(ur.df(residuals(lm(y~x-1)), lags=2))@teststat

  GHadf(y,x, kmax=12)$tstat
  summary(ur.df(residuals(lm(y~x-1)), lags=12))@teststat

}
##### PRINT METHOD
print.GregHan1996<-function(x,...){
  cat("#Gregory Hansen (1996) test of linear cointegration versus cointegration with breaks\n")
  cat("Model: 	", x$model,"\n")
  cat(paste("Number of right-side variables, m=", x$vars, COLLAPSE="\n"))
  cat("\n#ADF test \n")
  cat(" Lag choice method:",x$AdfLagChoice, "\t Value:", x$ADFlag, "\n")
  cat(" Test value:",x$adf, "\n")
  cat(" Critical values:\n")
  print(x$critical.values1)
  cat(" Breakpoint:", x$breakpointADF, "\n")

  cat("\n#Z-tau test \n")
  cat(" Test value:",x$Zt, "\n")
  cat(" Critical values:\n")
  print(x$critical.values1)
  cat(" Breakpoint:", x$breakpointZt, "\n")

  cat("\n#Z-alpha test \n")
  cat(" Test value:",x$Za, "\n")
  cat(" Critical values:\n")
  print(x$critical.values2)
  cat(" Breakpoint:", x$breakpointZa, "\n")
}

##### PLOT METHOD
plot.GregHan1996<-function(x,showTest=c("All", "ADF", "Zt","Za"),...)
{
  showTest<-match.arg(showTest)
  nplots<-if(showTest=="All") 3 else 1
  layout(matrix(1:nplots, nrow=nplots))
  forMode<-paste("for model ", x$model, sep="")
#ADF test plot
  if(showTest%in%c("All", "ADF")){
    temp1<-x$TestValues[, "temp1"]
    ts.plot(temp1, ylim=c(min(temp1,x$critical.values1[1])+0.1,max(temp1)+0.2), ylab="ADF", main="ADF statistic") 
    abline(h=x$critical.values1[1], col=2)
    abline(h=x$critical.values1[3], col=3)
    abline(h=x$critical.values1[4], col=4)
    legend("topleft", c("1%", "5%", "10%"), col=c(2,3,4), lty=c(1,1,1), ncol=3, box.lty=0)
  }	
  
#Zt test
  if(showTest%in%c("All", "Zt")){  
    temp4<-x$TestValues[, "temp4"]
    ts.plot(temp4, ylim=c(min(temp4,x$critical.values1[1])+0.1,max(temp1)+0.2), ylab="Zt", main="Zt statistic") 
    abline(h=x$critical.values1[1], col=2)
    abline(h=x$critical.values1[3], col=3)
    abline(h=x$critical.values1[4], col=4)
    legend("topleft", c("1%", "5%", "10%"), col=c(2,3,4), lty=c(1,1,1), ncol=3, box.lty=0)
  }

#Zt test
  if(showTest%in%c("All", "Za")){  
    temp3<-x$TestValues[, "temp3"]
    ts.plot(temp3, ylim=c(min(temp3,x$critical.values2[1])+0.1,max(temp1)+0.2), ylab="Za", main="Za statistic") 
    abline(h=x$critical.values2[1], col=2)
    abline(h=x$critical.values2[3], col=3)
    abline(h=x$critical.values2[4], col=4)
    legend("topleft", c("1%", "5%", "10%"), col=c(2,3,4), lty=c(1,1,1), ncol=3, box.lty=0)
  }
}




####Test of the function with internal data


#load data
if(!"MyresFin"%in%ls()) { load(file="~/Dropbox/Documents/stats/R/greg Han/MyresFin.rda")
  attach(MyresFin)
y<-md
x<-cbind(ynr,sr)
}



stTestC<-LRRbreak.test(y, x, model="C", AdfLagChoice="AIC", ADFlagmax=6)
sink(file=paste("~/Dropbox/Documents/stats/R/greg Han/GregHanResult", format(Sys.time(), "%a %b %d %H %Y"), ".txt"))
print(stTestC)
sink(file=NULL)

system(paste("diff", " '/home/mat/Dropbox/Documents/stats/R/greg Han/GregHanResult ", format(Sys.time(), "%a %b %d %H %Y"), " .txt'", " '/home/mat/Dropbox/Documents/stats/R/greg Han/GregHanResult Mon May  3 00:38:35 2010 .txt'",sep="" ), intern=TRUE)

plot(stTestC)
  
stTestCS<-LRRbreak.test(y, x, model="CS", AdfLagChoice="AIC", ADFlagmax=6)
print(stTestCS)
plot(stTestCS)
plot(stTestCS, showTest="ADF")

stTestCT<-LRRbreak.test(y, x, model="CT", AdfLagChoice="AIC", ADFlagmax=6)
print(stTestCT)
plot(stTestCT)

system.time(LRRbreak.test(y, x, model="CS", AdfLagChoice="AIC", ADFlagmax=6))
system.time(LRRbreak.test(y, x, model="CS", AdfLagChoice="AIC", ADFlagmax=6))
# Mon46m <- read.table (file="~/Documents/Ordi/MatLab/commandes/Gregory Hansen/mon46m.txt", header=FALSE, sep='\t', quote='"\'', dec='.', col.names = "Mon46qm")[,1]
# Mon46q<- read.table (file="~/Documents/Ordi/MatLab/commandes/Gregory Hansen/mon46q.txt", header=FALSE, sep='\t', quote='"\'', dec='.', col.names = "Mon46q")[,1]
# Mon47q<- read.table (file="~/Documents/Ordi/MatLab/commandes/Gregory Hansen/mon47q.txt", header=FALSE, sep='\t', quote='"\'', dec='.', col.names = "Mon47q")[,1]
