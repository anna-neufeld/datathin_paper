library("argparse")
library("tidyr")
source("GCS_splitfuns.R")
source("GCS_evalfuns.R")
library("dplyr")
library("mclust")
options(dplyr.summarise.inform = FALSE)

## -----------------------------------------
## load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--simname", default = "allsims",
                    help = "name of simulation")
parser$add_argument("--nreps", type = "double", default = 50,
                    help = "number of replicates for each set of params")
args <- parser$parse_args()
jobid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

nreps_per_combo <- args$nreps

## -----------------------------------------
## Binomial Simulation
## -----------------------------------------

filename <- paste("res/",args$simname, jobid, "_binomial.txt", sep="")

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

recon.bin.logit <- function(dat, lpsvd, N, k) {
  U <- lpsvd$u
  d <- lpsvd$d
  V <- lpsvd$v
  if (k == 1) {
    approx <- U[,1:k, drop='F']%*%(d[1:k])%*%t(V[,1:k, drop='F'])
  } else {
    approx <- U[,1:k]%*%diag(d[1:k])%*%t(V[,1:k])
  }
  logit.dat <- qlogis((dat+0.001)/(N+0.002))
  return(mean((logit.dat-approx)^2))
}

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

# Set the simulation parameters 
maxPC <- 20
maxf <- 5

n <- 250 
d <- 100

# Generate latent structure in the logit space
nPC <- 10 
sval <- 14:5 

allEps <- seq(0,1, length.out=51)[2:50]

for (i in 1:nreps_per_combo) {
  
  #### SET PARAMS
  Upop <- pracma::randortho(n)[,1:nPC]
  Dpop <- diag(sval)
  Vpop <- pracma::randortho(d)[,1:nPC]
  logitp <- Upop%*%Dpop%*%t(Vpop)
  N <- 100
  p <- plogis(logitp)
  
  
  # Generate binomial data
  X <- matrix(0, nrow=n, ncol=d)
  X[] <- rbinom(n*d, N, p)
  
  # Naive method
  p.hat <- (X+0.001)/(N+0.002)
  logitp.hat <- qlogis(p.hat)
  naive.svd <- svd(logitp.hat)
  
  for (j in 1:maxPC) {
    write(c("Naive", i, "RealMSE", j, recon.bin(X, naive.svd, N, j)), file = filename, append=TRUE, ncolumns = 5)
    write(c("Naive", i, "NLL", j, recon.nll(X, naive.svd, N, j)), file = filename, append=TRUE, ncolumns = 5)
  }
  
  # Data Thinning (0.5)
  eps <- 0.5
  sp <- ccsplit(X, family="binomial", epsilon=eps, arg=N)
  Ndt <- N*eps
  p.hat <- (sp$Xtr+0.001)/(Ndt+0.002)
  logitp.hat <- qlogis(p.hat)
  dt.svd <- svd(logitp.hat)
  
  for (j in 1:maxPC) {
    write(c("Data Thinning (0.5)", i, "RealMSE", j, recon.bin(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j)), file = filename, append=TRUE, ncolumns = 5)
    write(c("Data Thinning (0.5)", i, "NLL", j, recon.nll(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j)), file = filename, append=TRUE, ncolumns = 5)
  }
  
  # Data Thinning (0.8)
  eps <- 0.8
  sp <- ccsplit(X, family="binomial", epsilon=eps, arg=N)
  Ndt <- N*eps
  p.hat <- (sp$Xtr+0.001)/(Ndt+0.002)
  logitp.hat <- qlogis(p.hat)
  dt.svd <- svd(logitp.hat)
  
  for (j in 1:maxPC) {
    write(c("Data Thinning (0.8)", i, "RealMSE", j, recon.bin(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j)), file = filename, append=TRUE, ncolumns = 5)
    write(c("Data Thinning (0.8)", i, "NLL", j, recon.nll(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j)), file = filename, append=TRUE, ncolumns = 5)
  }
  
  # Multifold Data Thinning
  spcv <- cccv(X, family="binomial", nfolds=maxf, arg=N) %>%
    purrr::map(~mutate(as.data.frame(.x), idx = row_number())) %>% 
    bind_rows(.id="fold") %>%
    mutate(fold = as.numeric(fold))
  Nmft <- N*(maxf-1)/maxf
  
  mftreal <- matrix(0, nrow=maxf, ncol=maxPC)
  mftnll <- matrix(0, nrow=maxf, ncol=maxPC)
  for (f in 1:maxf) {
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
      mftreal[f,j] <- recon.bin(CVte, mft.svd, Nmft*1/(maxf-1), j)
      mftnll[f,j] <- recon.nll(CVte, mft.svd, Nmft*1/(maxf-1), j)
    }
  }
  
  for (j in 1:maxPC) {
    write(c("Multifold Thinning", i, "RealMSE", j, mean(mftreal[,j])), file = filename, append=TRUE, ncolumns = 5)
    write(c("Multifold Thinning", i, "NLL", j, mean(mftnll[,j])), file = filename, append=TRUE, ncolumns = 5)
  }
}

## -----------------------------------------
## Gamma Small Simulation
## -----------------------------------------

filename <- paste("res/",args$simname, jobid, "_gammasmall.txt", sep="")

n <- 100 
d <- 2
ncluster <- 4
true_clusters <- truth <- rep(c(1:ncluster), n)
theta0 <- 20
pars <- expand_grid(
  cluster = 1:ncluster,
  dim = 1:d
) %>% bind_cols(data.frame(shape = c(20,20,
                                     20,20,
                                     20,20,
                                     20,20), 
                           rate =  c(0.5,5,
                                     5,0.5,
                                     10,10,
                                     0.5,0.5)))
#### Constant for all simulations??? 
#### In future let's multiply this by things to change the signal strength. 
theta <- matrix(rep(pars$shape,n), nrow=n*ncluster, ncol=d, byrow=TRUE)
lambda <- matrix(rep(pars$rate,n), nrow=n*ncluster, ncol=d, byrow=TRUE)

maxk <- 10
maxf <- 5

for (trial in 1:nreps_per_combo) {
  #Generate Gamma data with the specified number of dimensions then split
  
  #### Ok, so theta and lambda were defined outside of this loop.
  #### Im ok with that. 
  Z <- matrix(0, nrow=n*ncluster, ncol=d)
  Z[] <- rgamma(n*d*ncluster, shape=theta, rate=lambda)
  
  sp <- ccsplit(Z, family="gamma", epsilon=0.5, arg=theta)
  sp2 <- ccsplit(Z, family="gamma", epsilon=0.8, arg=theta)
  spcv <- cccv(Z, family="gamma", nfolds=maxf, arg=theta) %>%
    purrr::map(~mutate(as.data.frame(.x), idx = row_number())) %>% 
    bind_rows(.id="fold") %>%
    mutate(fold = as.numeric(fold))
  
  trueCluster.holder <- kmeans(lambda, centers=ncluster)
  for (i in 1:ncluster) {
    trueCluster.holder$centers[i,] <- colMeans(sp$Xtr[trueCluster.holder$cluster==i,])
  }
  
  for (k in 1:maxk) {
    #First, consider the traditional within sum of squares approach
    sskmeans <- kmeans(Z, centers=k, nstart=10)
    rand <- mclust::adjustedRandIndex(sskmeans$cluster, true_clusters)
    write(c("Naive", rand, trial, "MSE", k, sskmeans$tot.withinss/(n*ncluster)), file = filename, append=TRUE, ncolumns = 6)
    write(c("Naive", rand, trial, "NLL", k, gammaNLL3(sskmeans, Z, Z, "naive", 0.5)), file = filename, append=TRUE, ncolumns = 6)
    
    #Then, consider the count splitting approach
    cskmeans <- kmeans(sp$Xtr, centers=k, nstart=10)
    rand <- mclust::adjustedRandIndex(cskmeans$cluster, true_clusters)
    write(c("Data Thinning (0.5)", rand, trial, "MSE", k, GCSMSE(cskmeans, sp$Xte, 0.5)), file = filename, append=TRUE, ncolumns = 6)
    write(c("Data Thinning (0.5)", rand, trial, "NLL", k, gammaNLL3(cskmeans, sp$Xtr, sp$Xte, "GCS", 0.5)), file = filename, append=TRUE, ncolumns = 6)
    
    #Then, consider the count splitting approach (eps = 0.8)
    cskmeans <- kmeans(sp2$Xtr, centers=k, nstart=10)
    rand <- mclust::adjustedRandIndex(cskmeans$cluster, true_clusters)
    write(c("Data Thinning (0.8)", rand, trial, "MSE", k, GCSMSE(cskmeans, sp2$Xte, 0.8)), file = filename, append=TRUE, ncolumns = 6)
    write(c("Data Thinning (0.8)", rand, trial, "NLL", k, gammaNLL3(cskmeans, sp2$Xtr, sp2$Xte, "GCS", 0.8)), file = filename, append=TRUE, ncolumns = 6)
    
    #Then, consider count splitting cross-validation strategy
    cvrand <- rep(0, maxf)
    cvmse <- rep(0, maxf)
    cvnll <- rep(0, maxf)
    for (f in 1:maxf) {
      CVtr <- spcv %>% filter(fold != f) %>% group_by(idx) %>% summarise(across(-fold, sum)) %>%
        select(-idx) %>% as.matrix
      CVte <- spcv %>% filter(fold == f) %>% select(-fold, -idx) %>% as.matrix
      cvkmeans <- kmeans(CVtr, centers=k, nstart=10)
      
      cvrand[f] <- mclust::adjustedRandIndex(cvkmeans$cluster, true_clusters)
      cvmse[f] <- GCSMSE(cvkmeans, CVte, (maxf-1)/maxf)
      cvnll[f] <- gammaNLL3(cvkmeans, CVtr, CVte, "GCS", (maxf-1)/maxf)
    }
    
    write(c("Multifold Thinning", mean(cvrand), trial, "MSE", k, mean(cvmse)), file = filename, append=TRUE, ncolumns = 6)
    write(c("Multifold Thinning", mean(cvrand), trial, "NLL", k, mean(cvnll)), file = filename, append=TRUE, ncolumns = 6)
  }
}

## -----------------------------------------
## Gamma Big Simulation
## -----------------------------------------

filename <- paste("res/",args$simname, jobid, "_gammabig.txt", sep="")

n <- 100
d <- 100

ncluster <- 10
true_clusters <- rep(c(1:ncluster), length.out=n*ncluster)
theta0 <- 2
theta <- matrix(rep(theta0, n*ncluster*d), nrow=n*ncluster, ncol=d, byrow=TRUE)
lambda <- matrix(rep(c(c(rep(0.1,20), rep(1,80)),
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
maxf <- 5

for (trial in 1:nreps_per_combo) {
  #Generate Gamma data with the specified number of dimensions then split
  
  #### Ok, so theta and lambda were defined outside of this loop.
  #### Im ok with that. 
  Z <- matrix(0, nrow=n*ncluster, ncol=d)
  Z[] <- rgamma(n*d*ncluster, shape=theta, rate=lambda)
  
  sp <- ccsplit(Z, family="gamma", epsilon=0.5, arg=theta)
  sp2 <- ccsplit(Z, family="gamma", epsilon=0.8, arg=theta)
  spcv <- cccv(Z, family="gamma", nfolds=maxf, arg=theta) %>%
    purrr::map(~mutate(as.data.frame(.x), idx = row_number())) %>% 
    bind_rows(.id="fold") %>%
    mutate(fold = as.numeric(fold))
  
  trueCluster.holder <- kmeans(lambda, centers=ncluster)
  for (i in 1:ncluster) {
    trueCluster.holder$centers[i,] <- colMeans(sp$Xtr[trueCluster.holder$cluster==i,])
  }
  
  for (k in 1:maxk) {
    #First, consider the traditional within sum of squares approach
    sskmeans <- kmeans(Z, centers=k, nstart=75)
    rand <- mclust::adjustedRandIndex(sskmeans$cluster, true_clusters)
    write(c("Naive", rand, trial, "MSE", k, sskmeans$tot.withinss/(n*ncluster)), file = filename, append=TRUE, ncolumns = 6)
    write(c("Naive", rand, trial, "NLL", k, gammaNLL3(sskmeans, Z, Z, "naive", 0.5)), file = filename, append=TRUE, ncolumns = 6)
    
    #Then, consider the count splitting approach
    cskmeans <- kmeans(sp$Xtr, centers=k, nstart=75)
    rand <- mclust::adjustedRandIndex(cskmeans$cluster, true_clusters)
    write(c("Data Thinning (0.5)", rand, trial, "MSE", k, GCSMSE(cskmeans, sp$Xte, 0.5)), file = filename, append=TRUE, ncolumns = 6)
    write(c("Data Thinning (0.5)", rand, trial, "NLL", k, gammaNLL3(cskmeans, sp$Xtr, sp$Xte, "GCS", 0.5)), file = filename, append=TRUE, ncolumns = 6)
    
    #Then, consider the count splitting approach (eps = 0.8)
    cskmeans <- kmeans(sp2$Xtr, centers=k, nstart=75)
    rand <- mclust::adjustedRandIndex(cskmeans$cluster, true_clusters)
    write(c("Data Thinning (0.8)", rand, trial, "MSE", k, GCSMSE(cskmeans, sp2$Xte, 0.8)), file = filename, append=TRUE, ncolumns = 6)
    write(c("Data Thinning (0.8)", rand, trial, "NLL", k, gammaNLL3(cskmeans, sp2$Xtr, sp2$Xte, "GCS", 0.8)), file = filename, append=TRUE, ncolumns = 6)
    
    #Then, consider count splitting cross-validation strategy
    cvrand <- rep(0, maxf)
    cvmse <- rep(0, maxf)
    cvnll <- rep(0, maxf)
    for (f in 1:maxf) {
      CVtr <- spcv %>% filter(fold != f) %>% group_by(idx) %>% summarise(across(-fold, sum)) %>%
        select(-idx) %>% as.matrix
      CVte <- spcv %>% filter(fold == f) %>% select(-fold, -idx) %>% as.matrix
      cvkmeans <- kmeans(CVtr, centers=k, nstart=75)
      
      cvrand[f] <- mclust::adjustedRandIndex(cvkmeans$cluster, true_clusters)
      cvmse[f] <- GCSMSE(cvkmeans, CVte, (maxf-1)/maxf)
      cvnll[f] <- gammaNLL3(cvkmeans, CVtr, CVte, "GCS", (maxf-1)/maxf)
    }
    
    write(c("Multifold Thinning", mean(cvrand), trial, "MSE", k, mean(cvmse)), file = filename, append=TRUE, ncolumns = 6)
    write(c("Multifold Thinning", mean(cvrand), trial, "NLL", k, mean(cvnll)), file = filename, append=TRUE, ncolumns = 6)
  }
}