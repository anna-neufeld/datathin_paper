nbSample <- function(x, bprime, eps) {
  p <- rbeta(length(x),eps*bprime,(1-eps)*bprime)
  return(rbinom(length(x),x,p))
}
normalSample <- function(x, sigmaprime, eps) {
  return(rnorm(length(x), eps*x, sqrt(eps*(1-eps)*sigmaprime)))
}

gammaSample <- function(x, alphaprime, eps) {
  Z <- rbeta(length(x), eps*alphaprime, (1-eps)*alphaprime)
  return(Z*x)
}


#### I was just verifying that the sampling functions are correct here. 

# n <- 100000
# normalmu <- 5
# normalsig <- 2
# X <- rnorm(n, normalmu, sqrt(normalsig))
# Xtrain <- normalSample(X, 3, 0.7)
# Xtest <- X-Xtrain
# cov(Xtrain, Xtest)
# cor(Xtrain, Xtest)
# 
# nbR <- 5
# nbP <- 0.7
# X <- rnbinom(n, nbR, nbP)
# Xtrain <- nbSample(X, 25, 0.3)
# Xtest <- X-Xtrain
# cor(Xtrain, Xtest)
# 
# alpha <- 5
# beta <- 7
# X <- rgamma(n, alpha, beta)
# Xtrain <- gammaSample(X, alpha, 0.3)
# Xtest <- X-Xtrain
# cor(Xtrain, Xtest)

exactCorExpression_nb <- function(ep,p, b, bprime) {
  varXtrain <- ep*(1-ep)*b*(1-p)/((bprime+1)*p)*(bprime+1/p+b*(1-p)/p)+ep^2*b*(1-p)/p^2
  varXtest <- ep*(1-ep)*b*(1-p)/((bprime+1)*p)*(bprime+1/p+b*(1-p)/p)+(1-ep)^2*b*(1-p)/p^2
  
  myCov <-   (ep*(1-ep)*b*(1-p)^2/p^2)*(1-(b+1)/(bprime+1))

  return(myCov/sqrt(varXtrain*varXtest))
}

exactCorExpression_normal <- function(ep,sigma, sigma.tilde) {
  varX <- sigma
  varXtrain <- ep^2*sigma+ep*(1-ep)*sigma.tilde
  varXtest <- (1-ep)^2*sigma+ep*(1-ep)*sigma.tilde
  #covar1 <- (varX - varXtrain - varXtest)/2
  covar2 <- ep*(1-ep)*(sigma-sigma.tilde)
  #cor1 <- covar1 / sqrt(varXtrain*varXtest)
  cor2 <- covar2 / sqrt(varXtrain*varXtest)
  return(cor2)
}

exactCorExpression_gamma <- function(ep, alpha, beta, alpha.tilde) {
  myVarXtrain <- ep*(1-ep)/(alpha.tilde+1)*(alpha/beta^2+alpha^2/beta^2)+ep^2*alpha/beta^2
  myVarXtest <- ep*(1-ep)/(alpha.tilde+1)*(alpha/beta^2+alpha^2/beta^2)+(1-ep)^2*alpha/beta^2
  myCov <- ep*(1-ep)*alpha/beta^2*(1-(alpha+1)/(alpha.tilde+1))

  return(myCov / sqrt(myVarXtrain*myVarXtest))
}


##### Negative Binomial Simulation 
n <- 100000
p <- 0.7
r <- 7
ep <- 0.44
tildeRs <- 10^(seq(-log10(100000), log10(1000000), length.out=50))
corsNB <- rep(NA, length(tildeRs ))
varXtrainNB <- rep(NA, length(tildeRs ))
varXtestNB <- rep(NA, length(tildeRs ))
c <- 1
X <- rnbinom(n, size=r, prob=p)

for (tildeR in tildeRs) {
  Xtrain <- sapply(X,function(u) nbSample(u, tildeR, ep))
  Xtest <- X-Xtrain
  corsNB[c] <-cor(Xtrain, Xtest)
  varXtrainNB[c] <- var(Xtrain)
  varXtestNB[c] <- var(Xtest)
  c <- c+1
}
trueCorsNB <- sapply(tildeRs, function(u) exactCorExpression_nb(ep,p,r,u))


##### Normal Simulation
mu <- 7
sigma <- 5
tildeSigmas <- 10^(seq(-log10(10), log10(1000), length.out=50))
corsNormal<- rep(NA, length(tildeSigmas))
varXtrainNormal <- rep(NA, length(tildeSigmas))
varXtestNormal <- rep(NA, length(tildeSigmas))
c <- 1
X <- rnorm(n, mu, sqrt(sigma))

for (tildeSigma in tildeSigmas) {
  Xtrain <- sapply(X,function(u) normalSample(u, tildeSigma, ep))
  Xtest <- X-Xtrain
  corsNormal[c] <-cor(Xtrain, Xtest)
  varXtrainNormal[c] <- var(Xtrain)
  varXtestNormal[c] <- var(Xtest)
  c <- c+1
}
trueCorsNormal <- sapply(tildeSigmas, function(u) exactCorExpression_normal(ep, sigma, u))

##### Gamma Simulation
alpha <- 7
beta <- 5
tildeAlphas <- 10^(seq(-log10(100), log10(1000), length.out=50))
corsGamma <- rep(NA, length(tildeAlphas))
varXtrainGamma <- rep(NA, length(tildeAlphas ))
varXtestGamma <- rep(NA, length(tildeAlphas ))
c <- 1
X <- rgamma(n,alpha, beta)

for (tildeAlpha in tildeAlphas) {
  Xtrain <- sapply(X,function(u) gammaSample(u, tildeAlpha, ep))
  Xtest <- X-Xtrain
  corsGamma[c] <-cor(Xtrain, Xtest)
  varXtrainGamma[c] <- var(Xtrain)
  varXtestGamma[c] <- var(Xtest)
  c <- c+1
}
trueCorsGamma <- sapply(tildeAlphas, function(u) exactCorExpression_gamma(ep, alpha, beta, u))

library(patchwork)
library(tidyverse)
pNB <- ggplot(data=NULL, aes(x=tildeRs, y=corsNB, col="Empirical"))+
  geom_point()+scale_x_log10()+
  geom_line(aes(x=tildeRs, y=trueCorsNB, col="Theoretical"))+labs(col="")+
  ylab("Correlation")+xlab(expression(tilde(r)~"(log scale)"))+
  ggtitle("Negative Binomial Distribution")

pNormal <- ggplot(data=NULL, aes(x=tildeSigmas, y=corsNormal, col="Empirical"))+
  geom_point()+scale_x_log10()+
  geom_line(aes(x=tildeSigmas, y=trueCorsNormal, col="Theoretical"))+labs(col="")+
  ylab("Correlation")+xlab(expression(tilde(sigma)^2~"(log scale)"))+
  ggtitle("Normal Distribution")

pGamma <- ggplot(data=NULL, aes(x=tildeAlphas, y=corsGamma, col="Empirical"))+
  geom_point()+scale_x_log10()+
  geom_line(aes(x=tildeAlphas, y=trueCorsGamma, col="Theoretical"))+labs(col="")+
  ylab("Correlation")+xlab(expression(tilde(alpha)~"(log scale)"))+
  ggtitle("Gamma Distribution")

pNormal+pNB+pGamma+plot_layout(guides="collect") & 
  theme_bw() & ylim(-1,1)
setwd("~/Dropbox/Generalized Count Splitting/")
ggsave("Figures/corPlot.png", width=12, height=4)
