library("tidyr")
library("dplyr")
library("tidyr")
library("patchwork")

# Helper functions

## Function to data thin from a Binomial(N, p) distribution 
binomsplit <- function(data, epsilon, pop) {
  X <- matrix(nrow=dim(data)[1], ncol=dim(data)[2])
  Y <- matrix(nrow=dim(data)[1], ncol=dim(data)[2])
  n <- length(data)
  
  if (is.null(pop)) {
    print("Binomial population missing.")
  }
  
  X[] <- rhyper(n, epsilon*pop, (1-epsilon)*pop, data)
  Y <- data - X
  
  return(list(Xtr = X, Xte = Y))
}

# Function to data thin from a Gamma(alpha, beta) distribution where alpha is the shape parameter
gammasplit <- function(data, epsilon, shape) {
  X <- matrix(nrow=dim(data)[1], ncol=dim(data)[2])
  Y <- matrix(nrow=dim(data)[1], ncol=dim(data)[2])
  n <- length(data)
  
  X[] <- rbeta(n, epsilon*shape, (1-epsilon)*shape)
  X <- X * data
  Y <- data - X
  
  return(list(Xtr = X, Xte = Y))
}

## Compute the binomial PCA rank-K reconstruction MSE loss
recon.bin <- function(dat, lpsvd, N, k) {
  U <- lpsvd$u
  d <- lpsvd$d
  V <- lpsvd$v
  if (k == 1) {
    approx <- U[,1:k, drop='F']%*%(d[1:k])%*%t(V[,1:k, drop='F'])
  } else {
    approx <- U[,1:k]%*%diag(d[1:k])%*%t(V[,1:k])
  }
  approx2 <- N * plogis(approx)
  return(mean((dat-approx2)^2))
}

## Compute the binomial PCA rank-K reconstruction NLL loss
recon.nll <- function(dat, lpsvd, N, k) {
  U <- lpsvd$u
  d <- lpsvd$d
  V <- lpsvd$v
  if (k == 1) {
    approx <- U[,1:k, drop='F']%*%(d[1:k])%*%t(V[,1:k, drop='F'])
  } else {
    approx <- U[,1:k]%*%diag(d[1:k])%*%t(V[,1:k])
  }
  return(sum(-dbinom(dat, N, plogis(approx), log=TRUE)))
}

## Compute the Gamma clustering MSE loss
GCSMSE <- function(km, testdat, epsilon) {
  cent <- km$centers * (1-epsilon)/epsilon
  labels <- km$cluster
  mse <- 0
  for (i in 1:nrow(testdat)) {
    mse <- mse + as.matrix(dist(rbind(testdat[i,], cent[labels[i],])))[1,2]^2
  }
  
  return(mse/(nrow(testdat)))
}

## Compute the Gamma clustering NLL loss
gammaNLL3 <- function(km, traindat, testdat, method, epsilon) {
  ests <- bind_cols(
    data.frame(cluster=km$cluster),
    as.data.frame(traindat)
  ) %>%
    pivot_longer(-cluster, names_to="dim") %>%
    group_by(cluster, dim) %>% 
    summarise(N=n(),
              scale=ifelse(N == 1, 1/value, mean(value*log(value))-mean(value)*mean(log(value))), 
              shape=ifelse(N == 1, 1, mean(value)/scale)) %>%
    mutate(scale = ifelse(N > 2, (N/(N-1))*scale, scale),
           shape = ifelse(N > 2, shape - (1/N)*(3*shape - (2/3)*(shape/(1+shape)) - (4/5)*(shape/((1+shape)^2))), shape)) %>% 
    mutate(rate = 1/scale) %>%
    select(-N, -scale)
  
  bind_cols(
    data.frame(cluster=km$cluster),
    as.data.frame(testdat)
  ) %>% 
    pivot_longer(-cluster, names_to="dim") %>% 
    left_join(ests, by=c("cluster", "dim")) %>% 
    mutate(nll = -dgamma(value, shape=shape*(1-epsilon)/epsilon, rate=rate, log=TRUE)) %>% 
    filter(!is.infinite(nll)) %>%
    summarise(nll = sum(nll)) %>%
    pull
}

# Set the number of simulations
nreps <- 2000




# BINOMIAL PCA SIMULATION

## Set the binomial simulation parameters 
maxPC <- 20
n <- 250 
d <- 100

## Generate latent structure in the logit space
nPC <- 10
sval <- 14:5

allEps <- seq(0,1, length.out=51)[2:50]
results_binom <- matrix(0, nrow=2*nreps*length(allEps)*maxPC, ncol=5)
colnames(results_binom ) <- c("eps", "sim", "measure", "k", "value")

## Run the simulations 
counter <- 1
for (i in 1:nreps) {
  
  # Set the binomial parameters
  Upop <- pracma::randortho(n)[,1:nPC]
  Dpop <- diag(sval)
  Vpop <- pracma::randortho(d)[,1:nPC]
  logitp <- Upop%*%Dpop%*%t(Vpop)
  N <- 100
  p <- plogis(logitp)
  
  
  # Generate binomial data
  X <- matrix(0, nrow=n, ncol=d)
  X[] <- rbinom(n*d, N, p)
  
  for (eps in allEps) {
    # Thin the data
    sp <- binomsplit(X, epsilon=eps, pop=N)
    
    # Compute the SVD in the logit space
    Ndt <- N*eps
    p.hat <- (sp$Xtr+0.001)/(Ndt+0.002)
    logitp.hat <- qlogis(p.hat)
    dt.svd <- svd(logitp.hat)
    
    for (j in 1:maxPC) {
      # Compute reconstruction errors and save results
      results_binom[counter,] <- c(eps, i, "RealMSE", j, recon.bin(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j))
      counter <- counter + 1
      results_binom[counter,] <- c(eps, i, "NLL", j, recon.nll(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j))
      counter <- counter + 1
    }
  }
}

## Summarise binomial simulation results
consRes_binom <- results_binom %>% 
  as.data.frame() %>%
  mutate(eps = as.numeric(eps)) %>%
  group_by(sim, eps, measure) %>%
  summarize(minK = which.min(value)) %>% 
  group_by(eps,measure) %>% 
  summarize(propCorrect = mean(minK==10))



# SMALL GAMMA CLUSTERING SIMULATION

# Set the small gamma clustering simulation parameters 
n <- 100
d <- 2
ncluster <- 4
true_clusters <- truth <- rep(c(1:ncluster), n)

## Gamma parameters for the "small d, small K" setting
pars <- data.frame(lambda = c(20,20,20,20,20,20,20,20), 
                   theta = c(0.5,5,5,0.5,10,10,0.5,0.5))
lambda <- matrix(rep(pars$lambda,n), nrow=ncluster*n, ncol=d, byrow=TRUE)
theta <- matrix(rep(pars$theta,n), nrow=ncluster*n, ncol=d, byrow=TRUE)

maxk <- 10
allEps <- seq(0,1, length.out=52)[2:51]

results_gammasmall <- matrix(0, nrow=2*nreps*length(allEps)*maxk, ncol=6)
colnames(results_gammasmall) <- c("eps", "sim", "method", "measure", "k", "value")

## Run the simulations 
counter <- 1
for (trial in 1:nreps) {
  #Generate Gamma data with the specified number of dimensions then split
  Z <- matrix(0, nrow=n*ncluster, ncol=d)
  Z[] <- rgamma(n*d*ncluster, shape=lambda, rate=theta)
  
  for (eps in allEps) {
    # Thin the data
    sp <- gammasplit(X, epsilon=0.5, shape=lambda)
    
    for (k in 1:maxk) {
      # Estimate k clusters, compute loss functions and save results
      cskmeans <- kmeans(sp$Xtr, centers=k, nstart=10)
      results_gammasmall[counter,] <- c(eps, trial, "GCP", "MSE", k, GCSMSE(cskmeans, sp$Xte, eps))
      counter <- counter + 1
      results_gammasmall[counter,] <- c(eps, trial, "GCP", "NLL", k, gammaNLL3(cskmeans, sp$Xtr, sp$Xte, "GCS", eps))
      counter <- counter + 1
    }
  }
}

## Summarise small gamma simulation results
consRes_gammasmall <- results_gammasmall %>% 
  as.data.frame() %>% 
  mutate(eps = as.numeric(eps)) %>%
  group_by(sim, eps, method, measure) %>%
  summarize(minK = which.min(value)) %>% 
  group_by(eps, method, measure) %>% 
  summarize(propCorrect = mean(minK==4))



# LARGE GAMMA CLUSTERING SIMULATION

# Set the large gamma clustering simulation parameters 
n <- 100
d <- 100
ncluster <- 10
true_clusters <- truth <- rep(c(1:ncluster), n)

## Gamma parameters for the "large d, large K" setting
lambda <- matrix(rep(2, n*ncluster*d), nrow=n*ncluster, ncol=d, byrow=TRUE)
theta <- matrix(rep(c(c(rep(0.1,20), rep(1,80)),
                       c(rep(1,10), rep(0.1,20), rep(1,70)),
                       c(rep(1,20), rep(0.1,20), rep(1,60)),
                       c(rep(1,30), rep(0.1,20), rep(1,50)),
                       c(rep(1,40), rep(0.1,20), rep(1,40)),
                       c(rep(1,50), rep(0.1,20), rep(1,30)),
                       c(rep(1,60), rep(0.1,20), rep(1,20)),
                       c(rep(1,70), rep(0.1,20), rep(1,10)),
                       c(rep(1,80), rep(0.1,20)),
                       c(rep(1,100))),
                     n), nrow=n*ncluster, ncol=d, byrow=TRUE)

maxk <- 15
allEps <- seq(0,1, length.out=52)[2:51]

results_gammalarge <- matrix(0, nrow=2*nreps*length(allEps)*maxk, ncol=6)
colnames(results_gammalarge) <- c("eps", "sim", "method", "measure", "k", "value")

## Run the simulations 
counter <- 1
for (trial in 1:nreps) {
  #Generate Gamma data with the specified number of dimensions then split
  Z <- matrix(0, nrow=n*ncluster, ncol=d)
  Z[] <- rgamma(n*d*ncluster, shape=lambda, rate=theta)
  
  for (eps in allEps) {
    # Thin the data
    sp <- gammasplit(X, epsilon=0.5, shape=lambda)
    
    for (k in 1:maxk) {
      # Estimate k clusters, compute loss functions and save results
      cskmeans <- kmeans(sp$Xtr, centers=k, nstart=75)
      results_gammalarge[counter,] <- c(eps, trial, "GCP", "MSE", k, GCSMSE(cskmeans, sp$Xte, eps))
      counter <- counter + 1
      results_gammalarge[counter,] <- c(eps, trial, "GCP", "NLL", k, gammaNLL3(cskmeans, sp$Xtr, sp$Xte, "GCS", eps))
      counter <- counter + 1
    }
  }
}

## Summarise large gamma simulation results
consRes_gammalarge <- results_gammalarge %>% 
  as.data.frame() %>%
  mutate(eps = as.numeric(eps)) %>%
  group_by(sim, eps, method, measure) %>%
  summarize(minK = which.min(value)) %>% 
  group_by(eps, method, measure) %>% 
  summarize(propCorrect = mean(minK==10)) 






# Create the Role of Epsilon NLL Plot (Figure 4)

pgamma_small <- ggplot(data=consRes_gammasmall %>% filter(measure=="NLL"), aes(x=eps, y=propCorrect))+
  geom_point()+
  geom_line()+
  ggtitle("Proportion of times we select correct K")+
  theme_bw()+
  xlab(expression(epsilon))+ylab("Proportion of datasets for which we selected K=4")+
  ylim(c(0,1)) +
  ggtitle("Gamma Clustering - Small")

pgamma_large <- ggplot(data=consRes_gammalarge %>% filter(measure=="NLL"), aes(x=eps, y=propCorrect))+
  geom_point()+
  geom_line()+
  ggtitle("Proportion of times we select correct K")+
  theme_bw()+
  xlab(expression(epsilon))+ylab("Proportion of datasets for which we selected K=10")+
  ylim(c(0,1)) +
  ggtitle("Gamma Clustering - Large")

pbinom <- ggplot(data=consRes_binom %>% filter(measure=="NLL"), aes(x=eps, y=propCorrect))+
  geom_point()+
  geom_line()+
  theme_bw()+
  xlab(expression(epsilon))+ylab("Proportion of datasets for which we selected K=10")+
  ylim(c(0,1)) +
  ggtitle("Binomial PCA")

pbinom+pgamma_small+pgamma_large



# Create the Role of Epsilon MSE Plot (Figure 8)

pgamma_small2 <- ggplot(data=consRes_gammasmall %>% filter(measure=="MSE"), aes(x=eps, y=propCorrect))+
  geom_point()+
  geom_line()+
  ggtitle("Proportion of times we select correct K")+
  theme_bw()+
  xlab(expression(epsilon))+ylab("Proportion of datasets for which we selected K=4")+
  ylim(c(0,1)) +
  ggtitle("Gamma Clustering - Small")

pgamma_large2 <- ggplot(data=consRes_gammalarge %>% filter(measure=="MSE"), aes(x=eps, y=propCorrect))+
  geom_point()+
  geom_line()+
  ggtitle("Proportion of times we select correct K")+
  theme_bw()+
  xlab(expression(epsilon))+ylab("Proportion of datasets for which we selected K=10")+
  ylim(c(0,1)) +
  ggtitle("Gamma Clustering - Large")

pbinom2 <- ggplot(data=consRes_binom %>% filter(measure=="RealMSE"), aes(x=eps, y=propCorrect))+
  geom_point()+
  geom_line()+
  theme_bw()+
  xlab(expression(epsilon))+ylab("Proportion of datasets for which we selected K=10")+
  ylim(c(0,1)) +
  ggtitle("Binomial PCA")


pbinom2+pgamma_small2+pgamma_large2


