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
  suppressMessages(
  {
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
  )
}
