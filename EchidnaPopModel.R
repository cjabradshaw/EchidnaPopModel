##########################################################################################################################################
## echidna (Tachyglossus aculeatus) demographic model
## 
## Corey Bradshaw
## corey.bradshaw@flinders.edu.au
## Flinders University, September 2021
##########################################################################################################################################

## functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

## source
source("matrixOperators.r")


##############################
## TACHYGLOSSUS (aculeatus) (TA)

# mass
TA.mass <- 4.0 # mean(c(3.8,3.4,3.7,(mean(c(3.9,7))))) # short-beaked echidna Tachyglossus aculeatus (Nicol & Andersen 2007)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
TA.rm.allom.pred <- 10^(0.6914 - (0.2622*log10(TA.mass*1000)))
TA.rm.pred <- 0.40 # rm/year = 0.40 (Schmidt-Nielsen K, Dawson T J, Crawford EC (1966) Temperature regulation in the echidna (Tachyglossus aculeatus) J Cell Physiol 67 : 63-72)
TA.lm.pred <- exp(TA.rm.pred)

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
TA.D.pred.up <- (10^(1.91 + (-1.02*log10(TA.mass))))/2 # divided by 2 for females only
TA.D.pred.lo <- (10^(-1.17 + (-0.76*log10(TA.mass))))/2 # divided by 2 for females only
TA.D.pred <- TA.D.pred.up # animals/km2 (checks out)

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
TA.age.max.allom <- round(10^(0.89 + (0.13*log10(TA.mass*1000))), 0)
TA.age.max.allom # underestimated (torpor, hibernation, low BMR)
TA.age.max <- 45 # (Nicol & Andersen 2007)

## age vector
TA.age.vec <- 0:TA.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
TA.F.allom.pred <- exp(2.719 - (0.211*log(TA.mass*1000)))/2 # divided by 2 for females
TA.F.allom.pred
TA.F.egg <- 1/2 # one egg/year; /2 for daughters
TA.F.pbreed <- 0.55 # females reproductively active (Nicol & Morrow 2012)
# 17 females produced 22 young over 7 years: 22/17/7 = 0.1849, or 0.1849/2 = 0.0924 (Rismiller & McKelvey 2000)
TA.F.pred <- TA.F.egg*TA.F.pbreed

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
TA.alpha.allom <- ceiling(exp(-1.34 + (0.214*log(TA.mass*1000))))
TA.alpha <- 3

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
TA.s.tran <- ln.a.s + b.s*log(TA.mass*1000) + log(1)
TA.s.ad.yr.allom <- exp(-exp(TA.s.tran))
TA.s.ad.yr.allom
TA.s.ad.yr <- mean(c(0.94,0.98))

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (1.05*TA.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 1.9 # rate of mortality decline (also known as bt)
a2 <- 1 - (TA.s.ad.yr) # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.13 # rate of mortality increase
longev <- TA.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
TA.Sx <- c(0.99*TA.s.ad.yr, 1 - qx)
plot(x, TA.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
TA.s.sd.vec <- 0.05*TA.Sx

## pre-breeding design with 0-1 survival in first row
TA.m.vec <- c(0, 0, 0, 0.5*TA.F.pred, 0.75*TA.F.pred, rep(TA.F.pred,(TA.age.max-4))) # 
TA.m.sd.vec <- 0.05*TA.m.vec
plot(TA.age.vec, TA.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
TA.m.dat <- data.frame(TA.age.vec, TA.m.vec)
param.init <- c(0.5, 4, -4)
TA.fit.logp <- nls(TA.m.vec ~ a / (1+(TA.age.vec/b)^c), 
                   data = TA.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
TA.fit.logp.summ <- summary(TA.fit.logp)
plot(TA.age.vec, TA.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
TA.age.vec.cont <- seq(0,max(TA.age.vec),1)
TA.pred.p.mm <- coef(TA.fit.logp)[1] / (1+(TA.age.vec.cont/coef(TA.fit.logp)[2])^coef(TA.fit.logp)[3])
#DN.pred.p.mm <- ifelse(DN.pred.p.m > 1, 1, DN.pred.p.m)
lines(TA.age.vec.cont, TA.pred.p.mm,lty=2,lwd=3,col="red")

## create matrix
TA.popmat <- matrix(data = 0, nrow=TA.age.max+1, ncol=TA.age.max+1)
diag(TA.popmat[2:(TA.age.max+1),]) <- TA.Sx[-1]
TA.popmat[TA.age.max+1,TA.age.max+1] <- TA.Sx[TA.age.max+1]
TA.popmat[1,] <- TA.pred.p.mm
colnames(TA.popmat) <- c(0:TA.age.max)
rownames(TA.popmat) <- c(0:TA.age.max)
TA.popmat.orig <- TA.popmat ## save original matrix

## matrix properties
max.lambda(TA.popmat.orig) ## 1-yr lambda
TA.lm.pred
max.r(TA.popmat.orig) # rate of population change, 1-yr
TA.ssd <- stable.stage.dist(TA.popmat.orig) ## stable stage distribution
plot(TA.age.vec, TA.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(TA.popmat.orig, TA.age.max) # reproductive value
TA.gen.l <- G.val(TA.popmat.orig, TA.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
TA.pop.found <- round(area*TA.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
TA.init.vec <- TA.ssd * TA.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*TA.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

TA.tot.F <- sum(TA.popmat.orig[1,])
TA.popmat <- TA.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
TA.n.mat <- matrix(0, nrow=TA.age.max+1,ncol=(t+1))
TA.n.mat[,1] <- TA.init.vec

## set up projection loop
for (i in 1:t) {
  TA.n.mat[,i+1] <- TA.popmat %*% TA.n.mat[,i]
}

TA.n.pred <- colSums(TA.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(TA.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
TA.K.max <- 1*TA.pop.found
TA.K.vec <- c(1, 0.25*TA.K.max, TA.K.max/2, 0.75*TA.K.max, TA.K.max) 
TA.red.vec <- c(1,0.993,0.97,0.93,0.8845)
plot(TA.K.vec, TA.red.vec,pch=19,type="b")
TA.Kred.dat <- data.frame(TA.K.vec, TA.red.vec)

# logistic power function a/(1+(x/b)^c)
TA.param.init <- c(1, TA.K.max, 2)
TA.fit.lp <- nls(TA.red.vec ~ a/(1+(TA.K.vec/b)^c), 
                 data = TA.Kred.dat,
                 algorithm = "port",
                 start = c(a = TA.param.init[1], b = TA.param.init[2], c = TA.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
TA.fit.lp.summ <- summary(TA.fit.lp)
plot(TA.K.vec, TA.red.vec, pch=19,xlab="N",ylab="reduction factor")
TA.K.vec.cont <- seq(1,2*TA.pop.found,1)
TA.pred.lp.fx <- coef(TA.fit.lp)[1]/(1+(TA.K.vec.cont/coef(TA.fit.lp)[2])^coef(TA.fit.lp)[3])
lines(TA.K.vec.cont, TA.pred.lp.fx, lty=3,lwd=3,col="red")

TA.a.lp <- coef(TA.fit.lp)[1]
TA.b.lp <- coef(TA.fit.lp)[2]
TA.c.lp <- coef(TA.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
TA.n.mat <- matrix(0, nrow=TA.age.max+1, ncol=(t+1))
TA.n.mat[,1] <- TA.init.vec
TA.popmat <- TA.popmat.orig

## set up projection loop
for (i in 1:t) {
  TA.totN.i <- sum(TA.n.mat[,i])
  TA.pred.red <- as.numeric(TA.a.lp/(1+(TA.totN.i/TA.b.lp)^TA.c.lp))
  diag(TA.popmat[2:(TA.age.max+1),]) <- (TA.Sx[-1])*TA.pred.red
  TA.popmat[TA.age.max+1,TA.age.max+1] <- 0
  TA.popmat[1,] <- TA.m.vec
  TA.n.mat[,i+1] <- TA.popmat %*% TA.n.mat[,i]
}

TA.n.pred <- colSums(TA.n.mat)
plot(yrs, TA.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=TA.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 500
itdiv <- iter/10

TA.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
TA.s.arr <- TA.m.arr <- array(data=NA, dim=c(t+1, TA.age.max+1, iter))

for (e in 1:iter) {
  TA.popmat <- TA.popmat.orig
  
  TA.n.mat <- matrix(0, nrow=TA.age.max+1,ncol=(t+1))
  TA.n.mat[,1] <- TA.init.vec
  
  for (i in 1:t) {
    # stochastic survival values
    TA.s.alpha <- estBetaParams(TA.Sx, TA.s.sd.vec^2)$alpha
    TA.s.beta <- estBetaParams(TA.Sx, TA.s.sd.vec^2)$beta
    TA.s.stoch <- rbeta(length(TA.s.alpha), TA.s.alpha, TA.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    TA.fert.stch <- rnorm(length(TA.popmat[,1]), TA.m.vec, TA.m.sd.vec)
    TA.m.arr[i,,e] <- ifelse(TA.fert.stch < 0, 0, TA.fert.stch)
    
    TA.totN.i <- sum(TA.n.mat[,i], na.rm=T)
    TA.pred.red <- TA.a.lp/(1+(TA.totN.i/TA.b.lp)^TA.c.lp)
    
    diag(TA.popmat[2:(TA.age.max+1),]) <- (TA.s.stoch[-1])*TA.pred.red
    TA.popmat[TA.age.max+1,TA.age.max+1] <- 0
    TA.popmat[1,] <- TA.m.arr[i,,e]
    TA.n.mat[,i+1] <- TA.popmat %*% TA.n.mat[,i]
    
    TA.s.arr[i,,e] <- TA.s.stoch * TA.pred.red
    
  } # end i loop
  
  TA.n.sums.mat[e,] <- ((as.vector(colSums(TA.n.mat))/TA.pop.found))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

TA.n.md <- apply(TA.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
TA.n.up <- apply(TA.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
TA.n.lo <- apply(TA.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

par(mfrow=c(1,3))
plot(yrs,TA.n.md,type="l", main = "", xlab="year", ylab="pN1", lwd=2, ylim=c(0.95*min(TA.n.lo),1.05*max(TA.n.up)))
lines(yrs,TA.n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,TA.n.up,lty=2,col="red",lwd=1.5)

TA.s.add <- TA.m.add  <- rep(0, TA.age.max+1)
for (m in 1:iter) {
  TA.s.add <- rbind(TA.s.add, TA.s.arr[ceiling(TA.gen.l):(t+1),,m])
  TA.m.add <- rbind(TA.m.add, TA.m.arr[ceiling(TA.gen.l):(t+1),,m])
}
TA.s.add <- TA.s.add[-1,]
TA.m.add <- TA.m.add[-1,]

TA.s.md <- apply(TA.s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
TA.s.up <- apply(TA.s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
TA.s.lo <- apply(TA.s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(TA.age.vec,TA.s.md,type="l", main = "", xlab="age", ylab="s", lwd=2, ylim=c(0.95*min(TA.s.lo),1.05*(max(TA.s.up))))
lines(TA.age.vec,TA.s.lo,lty=2,col="red",lwd=1.5)
lines(TA.age.vec,TA.s.up,lty=2,col="red",lwd=1.5)

TA.m.md <- apply(TA.m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
TA.m.up <- apply(TA.m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
TA.m.lo <- apply(TA.m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(TA.age.vec,TA.m.md,type="l", main = "", xlab="age", ylab="m", lwd=2, ylim=c(0.95*min(TA.m.lo),1.05*max(TA.m.up)))
lines(TA.age.vec,TA.m.lo,lty=2,col="red",lwd=1.5)
lines(TA.age.vec,TA.m.up,lty=2,col="red",lwd=1.5)
par(mfrow=c(1,1))

