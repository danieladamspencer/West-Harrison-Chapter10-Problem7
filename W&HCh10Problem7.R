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
m0.star <- c(-50,25,50,-25)
C0.star <- diag(c(625,225,625,225))

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
plot(sales,main="Sales of an Confectionary Product")
lines(ts(ft,frequency = 4,start=c(1975,1)),col='blue')
legend('bottomright',legend = c("Observed","Forecasted"),col=c('black','blue'),
       lty=c(1,1),bty='n')
