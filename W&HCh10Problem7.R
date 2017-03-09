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
