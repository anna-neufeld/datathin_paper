library("tidyr")
library("dplyr")
library("ggplot2")
library("patchwork")

source("simFunctions.R")


# Set the number of simulations
nreps <- 2000

## -----------------------------------------
## Binomial Simulation
## -----------------------------------------

# Set the simulation parameters 
maxPC <- 20
maxf <- 5

n <- 250 
d <- 100

# Generate latent structure in the logit space
nPC <- 10 
sval <- 14:5 

results_binom <- matrix(0, nrow=2*nreps*4*maxPC, ncol=5)
colnames(results_binom ) <- c("method", "sim", "measure", "k", "value")

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
  
  # Naive method that reuses data
  p.hat <- (X+0.001)/(N+0.002)
  logitp.hat <- qlogis(p.hat)
  naive.svd <- svd(logitp.hat)
  
  for (j in 1:maxPC) {
    # Compute reconstruction errors and save results
    results_binom[counter,] <- c("Naive", i, "MSE", j, recon.bin(X, naive.svd, N, j))
    counter <- counter + 1
    results_binom[counter,] <- c("Naive", i, "NLL", j, recon.nll(X, naive.svd, N, j))
    counter <- counter + 1
  }
  
  # Data Thinning (0.5)
  eps <- 0.5
  sp <- binomsplit(X, epsilon=eps, pop=N)
  Ndt <- N*eps
  p.hat <- (sp$Xtr+0.001)/(Ndt+0.002)
  logitp.hat <- qlogis(p.hat)
  dt.svd <- svd(logitp.hat)
  
  for (j in 1:maxPC) {
    # Compute reconstruction errors and save results
    results_binom[counter,] <- c("Data Thinning (0.5)", i, "MSE", j, recon.bin(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j))
    counter <- counter + 1
    results_binom[counter,] <- c("Data Thinning (0.5)", i, "NLL", j, recon.nll(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j))
    counter <- counter + 1
  }
  
  # Data Thinning (0.8)
  eps <- 0.8
  sp <- binomsplit(X, epsilon=eps, pop=N)
  Ndt <- N*eps
  p.hat <- (sp$Xtr+0.001)/(Ndt+0.002)
  logitp.hat <- qlogis(p.hat)
  dt.svd <- svd(logitp.hat)
  
  for (j in 1:maxPC) {
    # Compute reconstruction errors and save results
    results_binom[counter,] <- c("Data Thinning (0.8)", i, "MSE", j, recon.bin(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j))
    counter <- counter + 1
    results_binom[counter,] <- c("Data Thinning (0.8)", i, "NLL", j, recon.nll(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j))
    counter <- counter + 1
  }
  
  # Multifold Data Thinning
  spcv <- multibinom(X, maxf, arg=N) %>%
    purrr::map(~mutate(as.data.frame(.x), idx = row_number())) %>% 
    bind_rows(.id="fold") %>%
    mutate(fold = as.numeric(fold))
  Nmft <- N*(maxf-1)/maxf
  
  mftreal <- matrix(0, nrow=maxf, ncol=maxPC)
  mftnll <- matrix(0, nrow=maxf, ncol=maxPC)
  for (f in 1:maxf) {
    ## Reaggregate folds into train/test sets
    CVtr <- spcv %>% 
      filter(fold != f) %>% 
      group_by(idx) %>% 
      summarise(across(-fold, sum)) %>%
      select(-idx) %>% 
      as.matrix
    CVte <- spcv %>% 
      filter(fold == f) %>% 
      select(-fold, -idx) %>% 
      as.matrix
    
    p.hat <- (CVtr+0.001)/(Nmft+0.002)
    logitp.hat <- qlogis(p.hat)
    mft.svd <- svd(logitp.hat)
    
    for (j in 1:maxPC) {
      # Compute reconstruction errors and save results
      mftreal[f,j] <- recon.bin(CVte, mft.svd, Nmft*1/(maxf-1), j)
      mftnll[f,j] <- recon.nll(CVte, mft.svd, Nmft*1/(maxf-1), j)
    }
  }
  
  for (j in 1:maxPC) {
    # Compute reconstruction errors and save results
    results_binom[counter,] <- c("Multifold Thinning", i, "MSE", j, mean(mftreal[,j]))
    counter <- counter + 1
    results_binom[counter,] <- c("Multifold Thinning", i, "NLL", j, mean(mftnll[,j]))
    counter <- counter + 1
  }
}

## -----------------------------------------
## Gamma Small Simulation
## -----------------------------------------

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

results_gammasmall <- matrix(0, nrow=2*nreps*4*maxPC, ncol=5)
colnames(results_gammasmall) <- c("method", "sim", "measure", "k", "value")

## Run the simulations 
counter <- 1
for (trial in 1:nreps) {
  #Generate Gamma data with the specified number of dimensions then split
  Z <- matrix(0, nrow=n*ncluster, ncol=d)
  Z[] <- rgamma(n*d*ncluster, shape=lambda, rate=theta)
  
  sp <- gammasplit(Z, epsilon=0.5, shape=lambda)
  sp2 <- gammasplit(Z, epsilon=0.8, shape=lambda)
  spcv <- multigamma(Z, maxf, arg=lambda) %>%
    purrr::map(~mutate(as.data.frame(.x), idx = row_number())) %>% 
    bind_rows(.id="fold") %>%
    mutate(fold = as.numeric(fold))
  
  for (k in 1:maxk) {
    #First, consider the traditional within sum of squares approach
    sskmeans <- kmeans(Z, centers=k, nstart=10)
    results_gammasmall[counter,] <- c("Naive", trial, "MSE", k, sskmeans$tot.withinss/(n*ncluster))
    counter <- counter + 1
    results_gammasmall[counter,] <- c("Naive", trial, "NLL", k, gammaNLL3(sskmeans, Z, Z, "naive", 0.5))
    counter <- counter + 1
    
    #Then, consider the count splitting approach with eps = 0.5
    cskmeans <- kmeans(sp$Xtr, centers=k, nstart=10)
    results_gammasmall[counter,] <- c("Data Thinning (0.5)", trial, "MSE", k, GCSMSE(cskmeans, sp$Xte, 0.5))
    counter <- counter + 1
    results_gammasmall[counter,] <- c("Data Thinning (0.5)", trial, "NLL", k, gammaNLL3(cskmeans, sp$Xtr, sp$Xte, "GCS", 0.5))
    counter <- counter + 1
    
    #Then, consider the count splitting approach with eps = 0.8
    cskmeans <- kmeans(sp2$Xtr, centers=k, nstart=10)
    results_gammasmall[counter,] <- c("Data Thinning (0.8)", trial, "MSE", k, GCSMSE(cskmeans, sp2$Xte, 0.8))
    counter <- counter + 1
    results_gammasmall[counter,] <- c("Data Thinning (0.8)", trial, "NLL", k, gammaNLL3(cskmeans, sp2$Xtr, sp2$Xte, "GCS", 0.8))
    counter <- counter + 1
    
    #Then, consider count splitting cross-validation strategy
    cvmse <- rep(0, maxf)
    cvnll <- rep(0, maxf)
    for (f in 1:maxf) {
      CVtr <- spcv %>% filter(fold != f) %>% group_by(idx) %>% summarise(across(-fold, sum)) %>%
        select(-idx) %>% as.matrix
      CVte <- spcv %>% filter(fold == f) %>% select(-fold, -idx) %>% as.matrix
      cvkmeans <- kmeans(CVtr, centers=k, nstart=10)
      
      cvmse[f] <- GCSMSE(cvkmeans, CVte, (maxf-1)/maxf)
      cvnll[f] <- gammaNLL3(cvkmeans, CVtr, CVte, "GCS", (maxf-1)/maxf)
    }
    
    results_gammasmall[counter,] <- c("Multifold Thinning", trial, "MSE", k, mean(cvmse))
    counter <- counter + 1
    results_gammasmall[counter,] <- c("Multifold Thinning", trial, "NLL", k, mean(cvnll))
    counter <- counter + 1
  }
}

## -----------------------------------------
## Gamma Big Simulation
## -----------------------------------------

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

results_gammalarge <- matrix(0, nrow=2*nreps*4*maxPC, ncol=5)
colnames(results_gammalarge) <- c("method", "sim", "measure", "k", "value")

## Run the simulations 
counter <- 1
for (trial in 1:nreps) {
  #Generate Gamma data with the specified number of dimensions then split
  Z <- matrix(0, nrow=n*ncluster, ncol=d)
  Z[] <- rgamma(n*d*ncluster, shape=lambda, rate=theta)
  
  sp <- gammasplit(Z, epsilon=0.5, shape=lambda)
  sp2 <- gammasplit(Z, epsilon=0.8, shape=lambda)
  spcv <- multigamma(Z, maxf, arg=lambda) %>%
    purrr::map(~mutate(as.data.frame(.x), idx = row_number())) %>% 
    bind_rows(.id="fold") %>%
    mutate(fold = as.numeric(fold))
  
  for (k in 1:maxk) {
    #First, consider the traditional within sum of squares approach
    sskmeans <- kmeans(Z, centers=k, nstart=10)
    results_gammalarge[counter,] <- c("Naive", trial, "MSE", k, sskmeans$tot.withinss/(n*ncluster))
    counter <- counter + 1
    results_gammalarge[counter,] <- c("Naive", trial, "NLL", k, gammaNLL3(sskmeans, Z, Z, "naive", 0.5))
    counter <- counter + 1
    
    #Then, consider the count splitting approach with eps = 0.5
    cskmeans <- kmeans(sp$Xtr, centers=k, nstart=10)
    results_gammalarge[counter,] <- c("Data Thinning (0.5)", trial, "MSE", k, GCSMSE(cskmeans, sp$Xte, 0.5))
    counter <- counter + 1
    results_gammalarge[counter,] <- c("Data Thinning (0.5)", trial, "NLL", k, gammaNLL3(cskmeans, sp$Xtr, sp$Xte, "GCS", 0.5))
    counter <- counter + 1
    
    #Then, consider the count splitting approach with eps = 0.8
    cskmeans <- kmeans(sp2$Xtr, centers=k, nstart=10)
    results_gammalarge[counter,] <- c("Data Thinning (0.8)", trial, "MSE", k, GCSMSE(cskmeans, sp2$Xte, 0.8))
    counter <- counter + 1
    results_gammalarge[counter,] <- c("Data Thinning (0.8)", trial, "NLL", k, gammaNLL3(cskmeans, sp2$Xtr, sp2$Xte, "GCS", 0.8))
    counter <- counter + 1
    
    #Then, consider count splitting cross-validation strategy
    cvmse <- rep(0, maxf)
    cvnll <- rep(0, maxf)
    for (f in 1:maxf) {
      CVtr <- spcv %>% filter(fold != f) %>% group_by(idx) %>% summarise(across(-fold, sum)) %>%
        select(-idx) %>% as.matrix
      CVte <- spcv %>% filter(fold == f) %>% select(-fold, -idx) %>% as.matrix
      cvkmeans <- kmeans(CVtr, centers=k, nstart=10)
      
      cvmse[f] <- GCSMSE(cvkmeans, CVte, (maxf-1)/maxf)
      cvnll[f] <- gammaNLL3(cvkmeans, CVtr, CVte, "GCS", (maxf-1)/maxf)
    }
    
    results_gammalarge[counter,] <- c("Multifold Thinning", trial, "MSE", k, mean(cvmse))
    counter <- counter + 1
    results_gammalarge[counter,] <- c("Multifold Thinning", trial, "NLL", k, mean(cvnll))
    counter <- counter + 1
  }
}

save.image("res.RData")