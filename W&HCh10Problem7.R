##### Read in the Data #####
s <- c(157,227,240,191,
       157,232,254,198,
       169,234,241,167,
       163,227,252,185,
       179,261,264,196,
       179,248,256,193,
       186,260,270,210,
       171,227,241,170,
       140,218,208,193,
       184,235,245,209,
       206,260,264,227)
sales <- ts(s,frequency = 4,start = c(1975,1))
c <- c(10.6,8.5,6.7,4.1,
       1.9,0.4,1.1,1.9,
       0.2,2.9,4.1,-1.4,
       -4.0,-4.5,-5.3,-8.4,
       -12.8,-13.2,-10.1,-4.6,
       -1.1,-0.1,0.0,-2.5,
       -5.1,-6.4,-8.0,-6.5,
       -3.7,-1.3,6.1,16.5,
       22.9,23.9,18.0,8.3,
       2.9,0.7,-2.4,-7.0,
       -9.8,-10.6,-12.3,-13.2)
cost <- ts(c,frequency = 4,start=c(1975,1))

plot(sales,col='red')
plot(cost,col='blue')


##### Part A #####
#m0.star <- c(-50,25,50,-25)
#C0.star <- diag(c(625,225,625,225))

#DK CORRECTION offset by 1
m0.star <- c(-25,-50,25,50)
C0.star <- diag(c(225,625,225,625))

U <- c(crossprod(rep(1,4),C0.star)%*%rep(1,4))
A <- c(C0.star%*%rep(1,4)/U)

m0s <- c(m0.star - tcrossprod(A,rep(1,4))%*%m0.star)
C0s <- C0.star - tcrossprod(A,A)*U

library(xtable)
xCs <- xtable(C0s)
print(xCs,include.rownames = F)

m0 <- c(220,-1.5,m0s)
library(dlm)
C0 <- bdiag(matrix(c(225,0,0,0.49),2,2),C0s)
xC <- xtable(C0)
print(xC,include.rownames = F)

##### Part B #####
FF <- sapply(cost,function(z){c(1,z,1,rep(0,3))})
GG <- matrix(c(1,0,0,0,0,0,
               0,1,0,0,0,0,
               0,0,0,1,0,0,
               0,0,0,0,1,0,
               0,0,0,0,0,1,
               0,0,1,0,0,0),6,6,byrow = T)

a1 <- c(GG%*%m0)
GGCGG <- GG%*%C0%*%t(GG)

xGCG <- xtable(GGCGG)
print(xGCG,include.rownames = F)

##### Part C #####
TT <- length(sales)
dT <- dS <- 0.9
dR <- 0.98
S0 <- 100
n0 <- 12

at <- mt <- At <- matrix(NA,6,TT)
Rt <- Ct <- Wt <- array(NA,dim = c(TT,6,6))
ft <- qt <- St <- et <- nt <- numeric(TT)

at[,1] <- GG%*%m0
Wt[1,,] <- bdiag((1-dT)*C0[1,1]/dT,(1-dR)*C0[2,2]/dR,(1-dS)*C0[-(1:2),-(1:2)]/dS)
Rt[1,,] <- GG%*%C0%*%t(GG) + Wt[1,,]
qt[1] <- crossprod(FF[,1],Rt[1,,])%*%FF[,1] + S0
ft[1] <- crossprod(FF[,1],at[,1])
nt[1] <- n0 + 1
et[1] <- sales[1] - ft[1]
St[1] <- S0 + (S0/nt[1])*(et[1]^2 / qt[1] - 1)
At[,1] <- Rt[1,,]%*%FF[,1]/qt[1] 
Ct[1,,] <- (St[1]/S0)*(Rt[1,,] - tcrossprod(At[,1],At[,1])*qt[1])
mt[,1] <- at[,1] + At[,1]*et[1]

## Update equations ##
for(t in 2:TT){
  at[,t] <- GG%*%mt[,(t-1)]
  Wt[t,,] <- bdiag((1-dT)*Ct[(t-1),1,1]/dT,(1-dR)*Ct[(t-1),2,2]/dR,(1-dS)*Ct[(t-1),-(1:2),-(1:2)]/dS)
  Rt[t,,] <- GG%*%Ct[(t-1),,]%*%t(GG) + Wt[t,,]
  qt[t] <- crossprod(FF[,t],Rt[t,,])%*%FF[,t] + St[(t-1)]
  ft[t] <- crossprod(FF[,t],at[,t])
  nt[t] <- nt[(t-1)] + 1
  et[t] <- sales[t] - ft[t]
  St[t] <- S0 + (St[(t-1)]/nt[t])*(et[t]^2 / qt[t] - 1)
  At[,t] <- Rt[t,,]%*%FF[,t]/qt[t] 
  Ct[t,,] <- (St[t]/St[(t-1)])*(Rt[t,,] - tcrossprod(At[,t],At[,t])*qt[t])
  mt[,t] <- at[,t] + At[,t]*et[t]
}

## Plot showing the observed sales and the forecast values for each time
#DAN KIRSNER ADDITION, extract the quarterly components in an easier to understand fashion
q1=NULL
q2=NULL
q3=NULL
q4=NULL
for(m in 1:length(cost))
{
  q1=c(q1,mt[6-(m+3)%%4,m])
  q2=c(q2,mt[6-(m+0)%%4,m])
  q3=c(q3,mt[6-(m+1)%%4,m])
  q4=c(q4,mt[6-(m+2)%%4,m])
}


plot(sales,main="Sales of an Confectionary Product")
lines(ts(ft,frequency = 4,start=c(1975,1)),col='blue')
legend('bottomright',legend = c("Observed","Forecasted"),col=c('black','blue'),
       lty=c(1,1),bty='n')

## Plot showing the filtered estimates of the state vector
windows(21,14)
# pdf(file="~/GitHub/West-Harrison-Chapter10-Problem7/WriteUp/Filteredestimates.pdf",width = 21, height=14)
par(mfrow=c(2,3),cex.axis=2.5,cex.main=2)
plot(ts(mt[1,],frequency = 4,start=c(1975,1)),main="Trend Effect",ylab="",lwd=2)
plot(ts(mt[2,],frequency = 4,start=c(1975,1)),main="Regression Effect",ylab="",lwd=2)
plot(ts(q1,frequency = 4,start=c(1975,1)),main="First Quarter Effect",ylab="",lwd=2)
plot(ts(q2,frequency = 4,start=c(1975,1)),main="Second Quarter Effect",ylab="",lwd=2)
plot(ts(q3,frequency = 4,start=c(1975,1)),main="Third Quarter Effect",ylab="",lwd=2)
plot(ts(q4,frequency = 4,start=c(1975,1)),main="Fourth Quarter Effect",ylab="",lwd=2)
dev.off()

## Plot showing the forecast distribution for sales
fdist <- apply(cbind(ft,qt,nt),1,function(z){c(z[1]-qt(0.975,df=z[3])*sqrt(z[2]),z[1]+qt(0.975,df=z[3])*sqrt(z[2]),z[1])})
# pdf("~/GitHub/West-Harrison-Chapter10-Problem7/WriteUp/ForecastDist.pdf")
plot(sales,main="Forecast",ylim=c(min(c(fdist)),max(c(fdist))),lwd=3)
lines(ts(fdist[1,],frequency = 4,start = c(1975,1)),col='blue',lty=3,lwd=2)
lines(ts(fdist[2,],frequency = 4,start = c(1975,1)),col='blue',lty=3,lwd=2)
lines(ts(fdist[3,],frequency = 4,start = c(1975,1)),col='blue',lwd=3)
legend(1978,135,legend = c("Observed","Forecast","95% CI"),col=c("black",rep('blue',2)),lty=c(1,1,3),lwd=c(3,3,2))
dev.off()

##### Part D #####

## Plot showing the relative stability of the regression parameter over time
cdist <- apply(cbind(mt[2,],Ct[,2,2],nt),1,function(z){c(z[1] - qt(0.975,df=z[3])*sqrt(z[2]),z[1] + qt(0.975,df=z[3])*sqrt(z[2]))})
pdf("~/GitHub/West-Harrison-Chapter10-Problem7/WriteUp/CostandRegressionParm.pdf")
plot(cost,main="Cost and Regression Parameter",lwd=3,ylab='')
apply(cdist,1,function(z){lines(ts(z,frequency = 4,start = c(1975,1)),col='red',lty=3,lwd=2)})
lines(ts(mt[2,],frequency = 4,start = c(1975,1)),col='red',lwd=3)
legend('topleft',legend=c("Cost","Parameter Estimate","95% Parameter CI"),col=c('black',rep('red',2)),lwd=c(3,3,2),lty=c(1,1,3),bty='n')
dev.off()

#Refit model with non dynamic slope (by setting that discount factor to 1)
#then show predictions are the same
TT_non_dynamic <- length(sales)
dT_non_dynamic  <- dS_non_dynamic  <- 0.9
dR_non_dynamic  <- 1
S0_non_dynamic  <- 100
n0_non_dynamic  <- 12

at_non_dynamic  <- mt_non_dynamic  <- At_non_dynamic  <- matrix(NA,6,TT_non_dynamic )
Rt_non_dynamic  <- Ct_non_dynamic  <- Wt_non_dynamic  <- array(NA,dim = c(TT_non_dynamic ,6,6))
ft_non_dynamic  <- qt_non_dynamic  <- St_non_dynamic  <- et_non_dynamic  <- nt_non_dynamic  <- numeric(TT_non_dynamic )

at_non_dynamic [,1] <- GG%*%m0
Wt_non_dynamic[1,,] <- bdiag((1-dT_non_dynamic )*C0[1,1]/dT_non_dynamic ,(1-dR_non_dynamic )*C0 [2,2]/dR_non_dynamic ,(1-dS_non_dynamic )*C0[-(1:2),-(1:2)]/dS_non_dynamic )
Rt_non_dynamic[1,,] <- GG%*%C0%*%t(GG) + Wt_non_dynamic[1,,]
qt_non_dynamic[1] <- crossprod(FF[,1],Rt_non_dynamic[1,,])%*%FF[,1] + S0
ft_non_dynamic[1] <- crossprod(FF[,1],at[,1])
nt_non_dynamic[1] <- n0_non_dynamic  + 1
et_non_dynamic[1] <- sales[1] - ft_non_dynamic [1]
St_non_dynamic[1] <- S0_non_dynamic  + (S0_non_dynamic /nt_non_dynamic [1])*(et_non_dynamic [1]^2 / qt_non_dynamic [1] - 1)
At_non_dynamic[,1] <- Rt_non_dynamic [1,,]%*%FF[,1]/qt_non_dynamic [1] 
Ct_non_dynamic[1,,] <- (St_non_dynamic [1]/S0_non_dynamic )*(Rt_non_dynamic [1,,] - tcrossprod(At_non_dynamic [,1],At_non_dynamic [,1])*qt_non_dynamic [1])
mt_non_dynamic[,1] <- at[,1] + At[,1]*et[1]

## Update equations ##
for(t in 2:TT_non_dynamic){
  at_non_dynamic[,t] <- GG%*%mt_non_dynamic[,(t-1)]
  Wt_non_dynamic[t,,] <- bdiag((1-dT_non_dynamic)*Ct_non_dynamic[(t-1),1,1]/dT_non_dynamic,(1-dR_non_dynamic)*Ct_non_dynamic[(t-1),2,2]/dR_non_dynamic,(1-dS_non_dynamic)*Ct_non_dynamic[(t-1),-(1:2),-(1:2)]/dS_non_dynamic)
  Rt_non_dynamic[t,,] <- GG%*%Ct_non_dynamic[(t-1),,]%*%t(GG) + Wt_non_dynamic[t,,]
  qt_non_dynamic[t] <- crossprod(FF[,t],Rt_non_dynamic[t,,])%*%FF[,t] + St_non_dynamic[(t-1)]
  ft_non_dynamic[t] <- crossprod(FF[,t],at_non_dynamic[,t])
  nt_non_dynamic[t] <- nt_non_dynamic[(t-1)] + 1
  et_non_dynamic[t] <- sales[t] - ft_non_dynamic[t]
  St_non_dynamic[t] <- S0_non_dynamic + (St_non_dynamic[(t-1)]/nt_non_dynamic[t])*(et_non_dynamic[t]^2 / qt_non_dynamic[t] - 1)
  At_non_dynamic[,t] <- Rt_non_dynamic[t,,]%*%FF[,t]/qt_non_dynamic[t] 
  Ct_non_dynamic[t,,] <- (St_non_dynamic[t]/St_non_dynamic[(t-1)])*(Rt_non_dynamic[t,,] - tcrossprod(At_non_dynamic[,t],At_non_dynamic[,t])*qt_non_dynamic[t])
  mt_non_dynamic[,t] <- at_non_dynamic[,t] + At_non_dynamic[,t]*et_non_dynamic[t]
}

## Plot showing the observed sales and the forecast values for each time
plot(sales,main="Dynamic versus static slope in the DLM")
lines(ts(ft,frequency = 4,start=c(1975,1)),col='green')
lines(ts(ft_non_dynamic,frequency=4,start=c(1975,1)),col="blue")
legend('bottom',legend = c("Observed","Forecasted non dynamic slope","Forecasted dynamic slope"),col=c('black','green','blue'),
       lty=c(1,1),bty='n')
mtext(side=3,line=-1,text=paste0("correlation between the two forecasts=",round(cor(ft,ft_non_dynamic),4)))

##### Part E #####

#please check DK's math, but I the discount factor predictions can be written cutely (like in H&W page 199 for single component)
cost_future=c(8.4, 10.6, 7.2 ,13.0,-2.9 ,-0.7 ,-6.4, -7.0, -14.9 ,-15.9 ,-18.0 ,-22.3)

FF_future <- sapply(cost_future,function(z){c(1,z,1,rep(0,3))})
at0=mt[,44]
rt0=Ct[44,,]
Rt_future=array(0,dim=c(12,6,6))
at_future=matrix(0,6,12)
qt_future=rep(0,12)
ft_future=rep(0,12)
discounter=c(1/dT,1/dR,rep(1/dS,4))


Rt_future[1,,]=GG%*%rt0%*%t(GG)%*%diag(discounter)
at_future[,1]=GG%*%at0
qt_future[1]=t(FF_future[,1])%*%Rt_future[1,,]%*%FF_future[,1]+St[44]
ft_future[1]=t(FF_future[,1])%*%at_future[,1]
for(z in 2:12)
{
  Rt_future[z,,]=GG%*%rt0%*%t(GG)%*%(diag(discounter^z))
  at_future[,z]=GG%*%at_future[,z-1]
  qt_future[z]=t(FF_future[,z])%*%Rt_future[z,,]%*%FF_future[,z]+St[44]
  ft_future[z]=t(FF_future[,z])%*%at_future[,z]
}
lower_bound=ft_future-qt(.95, 44-6)*sqrt(qt_future)
upper_bound=ft_future+qt(.95, 44-6)*sqrt(qt_future)

plot(sales,main="Sales of an Confectionary Product",xlim=c(1975,1989),ylim=c(min(c(sales,lower_bound)),max(c(sales,upper_bound))))
lines(ts(lower_bound,frequency = 4,start=c(1986,1)),col='green')
lines(ts(ft_future,frequency=4,start=c(1986,1)),col="blue")
lines(ts(upper_bound,frequency = 4,start=c(1986,1)),col='green')
legend('topleft',legend = c("Observed","Forecasted Mean","Forecasted 90% Prediction Interval"),col=c('black','blue','green'),
       lty=c(1,1),bty='n')
