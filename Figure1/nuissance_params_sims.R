#### Datathin from a N(mu, sigma) distribution, assuming that sigma=sigmaprime. 
normalSample <- function(x, sigmaprime, eps) {
  return(rnorm(length(x), eps*x, sqrt(eps*(1-eps)*sigmaprime)))
}


### Datathin from a NB(r,p) distribution, assuming that r=rprime. 
nbSample <- function(x, rprime, eps) {
  p <- rbeta(length(x),eps*rprime,(1-eps)*rprime)
  return(rbinom(length(x),x,p))
}

#### Datathin from a Gamma(alpha, beta) distribution, assuming that alpha=alphaprime. 
gammaSample <- function(x, alphaprime, eps) {
  Z <- rbeta(length(x), eps*alphaprime, (1-eps)*alphaprime)
  return(Z*x)
}


### Computes the correlation between X1 and X2 when we datathin the N(mu, sigma) distribution,
#### assuming that sigma=sigmaprime. Suggested by Prop 1. 
exactCorExpression_normal <- function(ep,sigma, sigma.tilde) {
  varXtrain <- ep^2*sigma+ep*(1-ep)*sigma.tilde
  varXtest <- (1-ep)^2*sigma+ep*(1-ep)*sigma.tilde
  covar2 <- ep*(1-ep)*(sigma-sigma.tilde)
  cor2 <- covar2 / sqrt(varXtrain*varXtest)
  return(cor2)
}

### Computes the correlation between X1 and X2 when we datathin the NB(b,p) distribution
### assuming that b=bprime. Suggested by Prop 2. 
exactCorExpression_nb <- function(ep,p, b, bprime) {
  varXtrain <- ep*(1-ep)*b*(1-p)/((bprime+1)*p)*(bprime+1/p+b*(1-p)/p)+ep^2*b*(1-p)/p^2
  varXtest <- ep*(1-ep)*b*(1-p)/((bprime+1)*p)*(bprime+1/p+b*(1-p)/p)+(1-ep)^2*b*(1-p)/p^2
  
  myCov <-   (ep*(1-ep)*b*(1-p)^2/p^2)*(1-(b+1)/(bprime+1))

  return(myCov/sqrt(varXtrain*varXtest))
}



### Computes the correlation between X1 and X2 when we datathin the Gamma(alpha, beta) distribution,
#### assuming that alpha=alphaprime. Suggested by Prop 3. 
exactCorExpression_gamma <- function(ep, alpha, beta, alpha.tilde) {
  myVarXtrain <- ep*(1-ep)/(alpha.tilde+1)*(alpha/beta^2+alpha^2/beta^2)+ep^2*alpha/beta^2
  myVarXtest <- ep*(1-ep)/(alpha.tilde+1)*(alpha/beta^2+alpha^2/beta^2)+(1-ep)^2*alpha/beta^2
  myCov <- ep*(1-ep)*alpha/beta^2*(1-(alpha+1)/(alpha.tilde+1))
  return(myCov / sqrt(myVarXtrain*myVarXtest))
}

##### Normal Simulation. Creates the left panel of Figure 1. 
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
ep <- 0.44
trueCorsNormal <- sapply(tildeSigmas, function(u) exactCorExpression_normal(ep, sigma, u))


##### Negative Binomial Simulation. Creates the center panel of Figure 1. 
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


##### Gamma Simulation. Creates the right panel of Figure 1. 
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



#### Final plotting code.
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
  theme_bw() & ylim(-1,1) &
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14), plot.title=element_text(size=16), 
        legend.text=element_text(size=14))
setwd("~/Dropbox/Generalized Count Splitting/JMLR-resubmit-sep-2023-v2/Figures/")
ggsave("Figures/corPlot.png", width=13, height=4)
