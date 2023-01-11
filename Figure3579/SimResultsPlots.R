library("tidyr")
library("dplyr")
library("ggplot2")
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

# Function to multithin from a Binomial(N, p) distribution 
multibinom <- function(data, nfolds=5, eps=NULL, arg=NULL) {
  if (is.null(eps)) {
    eps <- rep(1/nfolds, nfolds)
  }
  
  arg2 <- arg
  output <- list()
  resdat <- data
  for (i in 1:(nfolds-1)) {
    epsfold <- eps[i]/sum(eps[i:length(eps)]) 
    temp <- binomsplit(resdat, epsfold, arg2)
    
    output <- append(output, list(i = temp$Xtr))
    resdat <- temp$Xte
    arg2 <- arg2 * (1-epsfold)
  }
  output <- append(output, list(nfolds = resdat))
  names(output) <- 1:nfolds
  return(output)
}

# Function to multithin from a Gamma(alpha, beta) distribution where alpha is the shape parameter
multigamma <- function(data, nfolds=5, eps=NULL, arg=NULL) {
  if (is.null(eps)) {
    eps <- rep(1/nfolds, nfolds)
  }
  
  arg2 <- arg
  output <- list()
  resdat <- data
  for (i in 1:(nfolds-1)) {
    epsfold <- eps[i]/sum(eps[i:length(eps)]) 
    temp <- gammasplit(resdat, epsfold, arg2)
    
    output <- append(output, list(i = temp$Xtr))
    resdat <- temp$Xte
    arg2 <- arg2 * (1-epsfold)
  }
  output <- append(output, list(nfolds = resdat))
  names(output) <- 1:nfolds
  return(output)
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


## -----------------------------------------
## Create the Plots
## -----------------------------------------

# Aggregate results from all three settings
results_all <-
  bind_rows(
    results_binom %>% as_tibble %>% mutate(setting = "Binomial PCA"),
    results_gammasmall %>% as_tibble %>% mutate(setting = "Gamma Clustering - Small"),
    results_gammalarge %>% as_tibble %>% mutate(setting = "Gamma Clustering - Large")
  ) %>%
  mutate(across(c("sim", "k", "value"), ~as.numeric(.)))

# Create the NLL as a function of k plot (Figure 3)
pvalNLL <- results_all %>%
  filter(measure == "NLL") %>%
  group_by(setting, method, k) %>%
  summarise(value = mean(value)) %>%
  mutate(value = (value-min(value))/(max(value)-min(value))) %>%
  ungroup %>%
  ggplot(aes(x=k, y=value, colour=method)) +
  geom_line() +
  geom_point(data=(results_all %>% 
                     filter(measure == "NLL") %>%
                     group_by(setting, method, k) %>% 
                     summarise(value = mean(value)) %>% 
                     mutate(value = (value-min(value))/(max(value)-min(value))) %>%
                     filter(value == min(value))
  ), 
  pch=1, size=3, stroke=2) +
  facet_wrap(~factor(setting, levels=c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large")), scales="free_x") +
  scale_colour_manual(name=NULL, 
                      values=c("Data Thinning (0.5)" = "#F8766D",
                               "Data Thinning (0.8)" = "#00BFC4",
                               "Multifold Thinning" = "#7CAE00",
                               "Naive" = "#C77CFF")) +
  geom_vline(data=data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
                             val = c(10, 4, 10)), aes(xintercept=val)) +
  xlab("K") + 
  ylab("") +
  theme(legend.position = "right") +
  theme_bw() +
  coord_cartesian(ylim=c(0,0.25))

pvalNLL

# Create the NLL histogram (Figure 5)
phistNLL <- results_all %>% 
  filter(measure == "NLL") %>%
  group_by(setting, sim, method) %>%
  filter(value == min(value)) %>%
  ungroup %>%
  filter(method %in% c("Data Thinning (0.8)", "Multifold Thinning")) %>%
  ggplot(aes(x=k, fill=method)) +
  geom_vline(data=data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
                             val = c(10, 4, 10)), aes(xintercept=val)) +
  geom_bar(aes(y = (..count..)/nreps), position=position_dodge2(preserve="single")) + 
  facet_wrap(~factor(setting, levels=c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large")), scales="free_x") +
  scale_fill_manual(name=NULL, 
                    values=c("Data Thinning (0.8)" = "#00BFC4",
                             "Multifold Thinning" = "#7CAE00")) +
  geom_blank(data=expand_grid(
    data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
               val = c(6, 1, 6)),
    data.frame(method = c("Data Thinning (0.8)", "Multifold Thinning"))),
    aes(x=val)) +
  geom_blank(data=expand_grid(
    data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
               val = c(16, 11, 16)),
    data.frame(method = c("Data Thinning (0.8)", "Multifold Thinning"))),
    aes(x=val)) +
  ylim(c(0,1)) +
  xlab("K") + 
  ylab("") +
  theme(legend.position = "right") +
  theme_bw()


phistNLL

# Create the MSE as a function of k plot (Figure 7)
pvalMSE <- results_all %>%
  filter(measure == "MSE") %>%
  group_by(setting, method, k) %>%
  summarise(value = mean(value)) %>%
  mutate(value = (value-min(value))/(max(value)-min(value))) %>%
  ungroup %>%
  ggplot(aes(x=k, y=value, colour=method)) +
  geom_line() +
  geom_point(data=(results_all %>% 
                     filter(measure == "MSE") %>%
                     group_by(setting, method, k) %>% 
                     summarise(value = mean(value)) %>% 
                     mutate(value = (value-min(value))/(max(value)-min(value))) %>%
                     filter(value == min(value))
  ), 
  pch=1, size=3, stroke=2) +
  facet_wrap(~factor(setting, levels=c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large")), scales="free_x") +
  scale_colour_manual(name=NULL, 
                      values=c("Data Thinning (0.5)" = "#F8766D",
                               "Data Thinning (0.8)" = "#00BFC4",
                               "Multifold Thinning" = "#7CAE00",
                               "Naive" = "#C77CFF")) +
  geom_vline(data=data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
                             val = c(10, 4, 10)), aes(xintercept=val)) +
  xlab("K") + 
  ylab("") +
  theme(legend.position = "right") +
  theme_bw() +
  coord_cartesian(ylim=c(0,0.25))

pvalMSE

# Create the MSE histogram (Figure 9)
phistMSE <- results_all %>% 
  filter(measure == "MSE") %>%
  group_by(setting, sim, method) %>%
  filter(value == min(value)) %>%
  ungroup %>%
  filter(method %in% c("Data Thinning (0.8)", "Multifold Thinning")) %>%
  ggplot(aes(x=k, fill=method)) +
  geom_vline(data=data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
                             val = c(10, 4, 10)), aes(xintercept=val)) +
  geom_bar(aes(y = (..count..)/nreps), position=position_dodge2(preserve="single")) + 
  facet_wrap(~factor(setting, levels=c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large")), scales="free_x") +
  scale_fill_manual(name=NULL, 
                    values=c("Data Thinning (0.8)" = "#00BFC4",
                             "Multifold Thinning" = "#7CAE00")) +
  geom_blank(data=expand_grid(
    data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
               val = c(6, 1, 6)),
    data.frame(method = c("Data Thinning (0.8)", "Multifold Thinning"))),
    aes(x=val)) +
  geom_blank(data=expand_grid(
    data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
               val = c(16, 11, 16)),
    data.frame(method = c("Data Thinning (0.8)", "Multifold Thinning"))),
    aes(x=val)) +
  ylim(c(0,1)) +
  xlab("K") + 
  ylab("") +
  theme(legend.position = "right") +
  theme_bw()

phistMSE
