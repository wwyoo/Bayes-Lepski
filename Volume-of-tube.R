#######################################################################################
#Author: William Weimin Yoo
#
#R code below constructs Bayesian credible bands using the volume of tube method.
#The codes will reproduce all simulation results found in the paper 
#It is to be used with the code for Bayes Lepski's method
#######################################################################################
#load libraries
library(splines)
library(MASS)
par(mai = c(0.5, 0.5, 0.1, 0.1))

##################################################################################
#true function
f0 <- function(x){
  2*x - x^3 + exp(-50*(x-1/2)^2) 
}

#----------------------------------------------------------------------------------------
#function to calculate derivatives of b-spline and form its design matrix, taken from bs()
bsprime <- function(x, derivs, df = NULL, knots = NULL, degree = 3, intercept = TRUE, Boundary.knots = range(x))
{
    nx <- names(x)
    x <- as.vector(x)
    Boundary.knots <- sort(Boundary.knots)
    ord <- 1L + (degree <- as.integer(degree))
    if(!is.null(df) && is.null(knots)){
      nIknots <- df - ord + (1L - intercept)
      knots <- if(nIknots > 0L) {
           knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
             2L)[-c(1L, nIknots + 2L)]
           stats::quantile(x, knots)
      }
    }
    Aknots <- sort(c(rep(Boundary.knots, ord), knots))

    basis <- splineDesign(Aknots, x, ord, derivs)
    n.col <- ncol(basis)
    dimnames(basis) <- list(nx, 1L:n.col)
    basis
}

#------------------------------------------------------------------------------------------
####################################################################################
#Posterior mean, variance and empirical Bayes
####################################################################################
#calculate sigma2hat for empirical Bayes
#Computes the empirical Bayes estimate for sigma
sigma2EB <- function(y, cmiddle, B, eta){
  ydiff <- y - B %*% eta
  sigma2 = (crossprod(ydiff) - crossprod(forwardsolve(cmiddle, crossprod(B, ydiff)))) / n
  return(sigma2)
}

#returns posterior mean and variance
pfmeanf <- function(y, cmiddle, B, b, Omegainv, eta){
    ans = forwardsolve(cmiddle, crossprod(B, y) + Omegainv %*% eta)
    pmean = crossprod(t(b), backsolve(t(cmiddle), ans))
    return(pmean)
}

pfvarf <- function(cmiddle, b){ 
    pSigma = crossprod(forwardsolve(cmiddle, t(b)))
    return(pSigma)
}

#----------------------------------------------------------------------------------------------
#computes quantiles for volume of tube method
betaf <- function(x, a, BB){
  bprimex <- bsprime(x, derivs = rep(1, length(x)), knots = a, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
  bx <- bs(x, knots = a, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
  A = crossprod(bprimex, bx) - crossprod(bx, bprimex)
  com = A %*% solve(BB, t(bx))
  num = crossprod(com, solve(BB, com))
  denom = crossprod(t(bx), solve(BB, t(bx)))
  beta = (num ^ (1 / 2)) / (denom ^ (3 / 2))
  return(beta)
}

winfty <- function(x){
  beta * exp(-x ^ 2 / 2)/(2 * pi) + 1 - pnorm(x) - gamma / 2
}

######################################################################################
#data generation
######################################################################################
n = 2000  #number of data points
nsam = 1000 #number of samples from posterior (for fixed width bands)
nrep = 1 #number of data replicates for coverage calculations 
x <- seq(0, 1, length = n)

q <- 4  #cubic splines

f0x <- f0(x)
sigma0 <- sqrt(0.1)  #true sigma
xi2 = 10  #prior variance
xiinv2 = 1/ xi2

gamma = 0.05
nnew <- 1000
xnew <- seq(from = 0, to = 1, length.out = nnew)

#Initialize
jmax <- 40
jmin <- 4
j <- jmax - 1

Fstore <- matrix(ncol = nnew)
maxPSDstore <- c()

f0newx <- f0(xnew)
coveragefix <- rep(0, nrep)
radiusfix <- rep(0, nrep)
coveragevar <- rep(0, nrep)
radiusvar <- rep(0, nrep)
coveragecorr <- rep(0, nrep)
radiuscorr <- rep(0, nrep)

#do some computations before the loop
#the J stays the same for all replications, and so fix J from the beginning and compute
#coverage by Monte Carlo replications

J <- jhat <- 31 #computed using the Bayes Lepski's method
N <- J - q #optimal in terms of mse

eta <- rep(1, J)
Omegainv <- xiinv2 * diag(rep(1, J))

a = seq(from = 0, to = 1, length.out = N + 2)[-c(1, N + 2)]
B <- bs(x, knots = a, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
BB <- crossprod(B)

ptm <- proc.time()
#compute quantiles
#volume of tube
 beta <- mean(sapply(xnew, betaf, a = a, BB = BB))  #Monte Carlo integration
 wgamma <- uniroot(winfty, c(0, 10))[[1]]
vot.time <- proc.time() - ptm

ptm <- proc.time()
#direct sampling
 V = BB + Omegainv
 cmiddle <- t(chol(V))

 bnewx <- predict(B, xnew)
 pfvar <- bnewx %*% solve(V, t(bnewx))
 sd1 <- matrix(sqrt(diag(pfvar)), nnew, nnew, byrow = FALSE)
 sd2 <- matrix(sqrt(diag(pfvar)), nnew, nnew, byrow = TRUE)

 pfcorr <- pfvar / (sd1 * sd2)

 pfhat <- mvrnorm(n = 1000, mu = rep(0, nnew), Sigma = pfcorr)
 pzmax = apply(abs(pfhat), 1, max)
 pzmaxquan = quantile(pzmax, 1 - gamma)
direct.time <- proc.time() - ptm

ptm <- proc.time()
#fix ball
 pfhat <- mvrnorm(n = 1000, mu = rep(0, nnew), Sigma = pfvar)
 prmax = apply(abs(pfhat), 1, max)
 prmaxquan = quantile(prmax, 1 - gamma) 
fix.time <- proc.time() - ptm

set.seed(100)
for(i in 1:nrep){
 y = f0x + sigma0 * rnorm(n)  #generate some data
 pfmean = pfmeanf(y, cmiddle = cmiddle, B = B, Omegainv = Omegainv, eta = eta, b = bnewx)
 psigma2 <- as.numeric(sigma2EB(y, cmiddle = cmiddle, B = B, eta = eta))

################################################################################
#Volume of Tube method
################################################################################
 pradius <- wgamma * sqrt(psigma2 * diag(pfvar))

 cbandup <- pfmean + pradius
 cbandlower <- pfmean - pradius

 radiusvar[i] <- mean(pradius)

 countvar <- mean(cbandlower <= f0newx & f0newx <= cbandup) #band actual coverage

 if(countvar == 1){
  coveragevar[i] <- 1
 }
 if(countvar != 1){
  coveragevar[i] <- 0
 }

 plot(xnew, f0newx, type = "l", lty = 1, ylim = c(-1, 3), lwd = 5)
 lines(x, y, type = "p")
 lines(xnew, pfmean, type = "l", lty = 5, lwd = 5, col=4)
 lines(xnew, cbandup, type = "l", lty = 5, lwd = 5, col=2)
 lines(xnew, cbandlower, type = "l", lty = 5, lwd = 5, col=2)
 print(i)
}

############################################################################
#Direct sampling
############################################################################
 fcorrradius <- pzmaxquan * sqrt(psigma2 * diag(pfvar))

 radiuscorr[i] <- mean(fcorrradius)

 pfcorrbandup <- pfmean + fcorrradius
 pfcorrbandlower <- pfmean - fcorrradius

 countcorr <- mean(pfcorrbandlower <= f0newx & f0newx <= pfcorrbandup)

 if(countcorr == 1){
  coveragecorr[i] <- 1
 }
 if(countcorr != 1){
  coveragecorr[i] <- 0
 }

 plot(xnew, f0newx, type = "l", lty = 1, ylim = c(-1, 3), lwd = 5)
 lines(x, y, type = "p")
 lines(xnew, pfmean, type = "l", lty = 3, lwd = 5, col=4)
 lines(xnew, pfcorrbandup, type = "l", lty = 5, lwd = 5, col=2)
 lines(xnew, pfcorrbandlower, type = "l", lty = 5, lwd = 5, col=2)

 print(i)
}

#############################################################################
#Sup-norm ball around posterior mean with fixed (not inflated) radius
#############################################################################
# fradius <- sqrt(psigma2) * prmaxquan 

# radiusfix[i] <- fradius

# pfbandup <- pfmean + fradius
# pfbandlower <- pfmean - fradius

# countfix <- mean(pfbandlower <= f0newx & f0newx <= pfbandup)  

# if(countfix == 1){
#  coveragefix[i] <- 1
# }
# if(countfix != 1){
#  coveragefix[i] <- 0
# }

# plot(xnew, f0newx, type = "l", lty = 1, ylim = c(-1, 3), lwd = 5)
# lines(x, y, type = "p")
# lines(xnew, pfmean, type = "l", lty = 3, lwd = 5, col=4)
# lines(xnew, pfbandup, type = "l", lty = 5, lwd = 5, col=2)
# lines(xnew, pfbandlower, type = "l", lty = 5, lwd = 5,col=2)
# print(i)
#}

#mean(coveragevar)
#mean(coveragefix)

