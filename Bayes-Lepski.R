##########################################################################################
#Author: William Weimin Yoo
#
#The R code below chooses J the number of basis functions using the Bayes Lepski's method
#############################################################################################

#load libraries
library(splines)
library(MASS)

f0 <- function(x){
  2*x - x^3 + exp(-50*(x-1/2)^2) 
}

#Computes the empirical Bayes estimate for sigma
sigma2EB <- function(y, cmiddle, B, eta){
  ydiff <- y - B %*% eta
  sigma2 = (crossprod(ydiff) - crossprod(forwardsolve(cmiddle, crossprod(B, ydiff)))) / n
  return(sigma2)
}

#Computes posterior mean and max standard devation at level j
pmeanf <- function(j, y){
 N <- j - q
 eta <- rep(0, j)
 Omegainv <- xiinv2 * diag(rep(1, j))
 a <- seq(from = 0, to = 1, length.out = N + 2)[-c(1, N + 2)]
 B <- bs(x, knots = a, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
 BB <- crossprod(B)
 V <- BB + Omegainv
 cmiddle <- t(chol(V))
 b <- predict(B, xnew)

 ans = forwardsolve(cmiddle, crossprod(B, y) + Omegainv %*% eta)
 pmean = crossprod(t(b), backsolve(t(cmiddle), ans))

 pvar <- b %*% solve(V, t(b)) 
 sigmaEB <- sqrt(as.numeric(sigma2EB(y = y, cmiddle = cmiddle, B = B, eta = eta)))

 maxpsd <- sigmaEB * sqrt(min(diag(pvar)) * log(j))
 return(list("pmean" = pmean, "maxpsd" = maxpsd))
}

#Is the sup norm difference of posterior means <= max posterior standard deviation? 
compare <- function(j, F, maxPSD, y){
 rule <- 0
 pmeanj <- pmeanf(j = j, y = y)
 fj <- as.numeric(pmeanj$pmean)
 maxpsdj <- pmeanj$maxpsd

 rule <- sum( apply(  abs(  t(fj - t(F) )   ), 1,  max) > maxPSD )
 return(list("rule" = rule, "pmeanj" = fj, "maxpsdj" = maxpsdj))
}

#previous version with for loop
#compare <- function(j, F, maxPSD, y){
# rule <- 0
# pmeanj <- pmeanf(j = j, y = y)
# fj <- pmeanj$pmean
# maxpsdj <- pmeanj$maxpsd

# for(i in 1:nrow(F)){
#  rule <- rule + ( max(abs(fj - F[i, ])) > maxPSD[i] )
# }
# return(list("rule" = rule, "pmeanj" = fj, "maxpsdj" = maxpsdj))
#}

##################################################################################
#Generate data
set.seed(100)
n = 2000  #number of data points
x <- seq(0, 1, length = n)

q <- 4  #cubic splines

f0x = f0(x)
sigma0 <- sqrt(0.1)  #true sigma

xi2 = 10
xiinv2 = (1 / xi2)
nnew <- 1000
xnew <- seq(from = 0, to = 1, length.out = nnew)

#Initialize
jmax <- 40
jmin <- 4
j <- jmax - 1

Fstore <- matrix(ncol = nnew)
maxPSDstore <- c()

nrep <- 1000
jdis <- rep(0, nrep)

for(i in 1:nrep){
 y = f0x + sigma0 * rnorm(n)  #generate some data

 #Start from jmax
 fmax <- pmeanf(j = jmax, y = y)
 Fstore[1, ] <- fmax$pmean
 maxPSDstore[1] <- fmax$maxpsd

 while(j >= jmin){
  reduce <- compare(j = j, F = Fstore, maxPSD = maxPSDstore, y = y)
  if(reduce$rule > 0){
   jhat <- j + 1
   break
  }
  else{
   Fstore <- rbind(Fstore, c(reduce$pmeanj))
   maxPSDstore <- c(maxPSDstore, reduce$maxpsdj)
   j <- j - 1
  }
 }
#jhat
 jdis[i] <- jhat
 print(i)
}
jdis
