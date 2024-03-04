library("argparse")
library("tidyr")
source("GCS_splitfuns.R")
source("GCS_evalfuns.R")
library("dplyr")
library("tidyr")
library("mclust")

## -----------------------------------------
## load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--simname", default = "test",
                    help = "name of simulation")
parser$add_argument("--nreps", type = "double", default = 50,
                    help = "number of replicates for each set of params")
args <- parser$parse_args()
jobid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

nreps_per_combo <- args$nreps

filename <- paste("res/",args$simname, jobid, ".txt", sep="")

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
nPC <- 10 #11
sval <- 14:5 #15:5

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
  
  for (eps in allEps) {
    sp <- ccsplit(X, family="binomial", epsilon=eps, arg=N)
    Ndt <- N*eps
    p.hat <- (sp$Xtr+0.001)/(Ndt+0.002)
    logitp.hat <- qlogis(p.hat)
    dt.svd <- svd(logitp.hat)
  
    for (j in 1:maxPC) {
      write(c(eps, i, "RealMSE", j, recon.bin(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j)),file = filename, append=TRUE, ncolumns = 5)
      write(c(eps, i, "LogitMSE", j, recon.bin.logit(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j)) ,file = filename, append=TRUE, ncolumns = 5)
      write(c(eps, i, "NLL", j, recon.nll(sp$Xte, dt.svd, Ndt*(1-eps)/eps, j)) ,file = filename, append=TRUE, ncolumns = 5)
    }
  }
}
    
    
    
    
  
    
    
    
    
    


  