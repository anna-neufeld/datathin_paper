library("argparse")
library("tidyr")
source("GCS_splitfuns.R")
source("GCS_evalfuns.R")
library("dplyr")
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
#Z <- matrix(0, nrow=ncluster*n, ncol=d)
#Z[] <- rgamma(n*d*ncluster, shape=theta, rate=lambda)


maxk <- 15


options(dplyr.summarise.inform = FALSE)
allEps <- seq(0,1, length.out=52)[2:51]

for (trial in 1:nreps_per_combo) {
  #Generate Gamma data with the specified number of dimensions then split
  
  #### Ok, so theta and lambda were defined outside of this loop.
  #### Im ok with that. 
  Z <- matrix(0, nrow=n*ncluster, ncol=d)
  Z[] <- rgamma(n*d*ncluster, shape=theta, rate=lambda)
  
  
  for (eps in allEps) {
    sp <- ccsplit(Z, family="gamma", epsilon=eps, arg=theta)
    
    trueCluster.holder <- kmeans(lambda, centers=ncluster)
    for (i in 1:ncluster) {
      trueCluster.holder$centers[i,] <- colMeans(sp$Xtr[trueCluster.holder$cluster==i,])
    }
    
    #### If we want this we need
    lucyLossMSE <- GCSMSE(trueCluster.holder, sp$Xte, eps)
    lucyLossNLL <- gammaNLL3(trueCluster.holder, sp$Xtr, sp$Xte, "GCS", eps)
    for (k in 1:maxk) {
      #Then, consider the count splitting approach
      cskmeans <- kmeans(sp$Xtr, centers=k, nstart=75)
      rand <- mclust::adjustedRandIndex(cskmeans$cluster, true_clusters)
      write(c(eps, rand, lucyLossMSE, trial, "GCP", "MSE", k, GCSMSE(cskmeans, sp$Xte, eps)), file = filename, append=TRUE, ncolumns = 8)
      write(c(eps, rand, lucyLossNLL, trial, "GCP", "NLL", k, gammaNLL3(cskmeans, sp$Xtr, sp$Xte, "GCS", eps)), file = filename, append=TRUE,ncolumns = 8)
    }
  }
}


