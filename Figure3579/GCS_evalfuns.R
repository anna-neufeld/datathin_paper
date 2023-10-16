GCSMSE <- function(km, testdat, epsilon) {
  cent <- km$centers * (1-epsilon)/epsilon
  labels <- km$cluster
  mse <- 0
  for (i in 1:nrow(testdat)) {
    mse <- mse + as.matrix(dist(rbind(testdat[i,], cent[labels[i],])))[1,2]^2
  }
  
  return(mse/(nrow(testdat)))
}

SSMSE <- function(km, testdat, testlabs) {
  cent <- km$centers
  mse <- 0
  for (i in 1:nrow(testdat)) {
    mse <- mse + as.matrix(dist(rbind(testdat[i,], cent[testlabs[i],])))[1,2]^2
  }
  
  return(mse/nrow(testdat))
}

gammaMLE <- function(dat) {
  res <- c(0,0)
  if (length(dat) > 1) {
    res <- tryCatch(MASS::fitdistr(dat, "gamma")$estimate,
                    error = function(e) {return(MASS::fitdistr(dat, "gamma", lower=c(0,0), start=list(shape=10*mean(dat), rate=10))$estimate)})
  } else {
    res <- c(1, MASS::fitdistr(dat, "exponential")$estimate) #???
  }
  return(tibble(shape=res[1], rate=res[2]))
}

gammaNLL <- function(km, traindat, testdat, method, epsilon) {
  MLEs <- bind_cols(
    data.frame(cluster=km$cluster),
    as.data.frame(traindat)
  ) %>%
    pivot_longer(-cluster, names_to="dim") %>%
    group_by(cluster, dim) %>% 
    summarise(suppressWarnings(gammaMLE(value)))
  
  if (method != "SS") {
    bind_cols(
      data.frame(cluster=km$cluster),
      as.data.frame(testdat)
    ) %>% 
      pivot_longer(-cluster, names_to="dim") %>% 
      left_join(MLEs, by=c("cluster", "dim")) %>% 
      mutate(nll = -dgamma(value, shape=shape* (1-epsilon)/epsilon, rate=rate, log=TRUE)) %>% 
      filter(!is.infinite(nll)) %>%
      summarise(nll = sum(nll)) %>%
      pull
  } else {
    testdat %>% 
      pivot_longer(-cluster, names_to="dim") %>% 
      left_join(MLEs, by=c("cluster", "dim")) %>% 
      mutate(nll = -dgamma(value, shape=shape* (1-epsilon)/epsilon, rate=rate, log=TRUE)) %>% 
      summarise(nll = sum(nll)) %>%
      pull
  }
}

gammaNLL2 <- function(km, traindat, testdat, method, epsilon) {
  MLEs <- bind_cols(
    data.frame(cluster=km$cluster),
    as.data.frame(traindat)
  ) %>%
    pivot_longer(-cluster, names_to="dim") %>%
    group_by(cluster, dim) %>% 
    summarise(N=n(),
              scale=ifelse(N == 1, 1/value, mean(value*log(value))-mean(value)*mean(log(value))), 
              shape=ifelse(N == 1, 1, mean(value)/scale)) %>%
    mutate(rate = 1/scale) %>%
    select(-N, -scale)
  
  bind_cols(
    data.frame(cluster=km$cluster),
    as.data.frame(testdat)
  ) %>% 
    pivot_longer(-cluster, names_to="dim") %>% 
    left_join(MLEs, by=c("cluster", "dim")) %>% 
    mutate(nll = -dgamma(value, shape=shape*(1-epsilon)/epsilon, rate=rate, log=TRUE)) %>% 
    filter(!is.infinite(nll)) %>%
    summarise(nll = sum(nll)) %>%
    pull
}

gammaNLL3 <- function(km, traindat, testdat, method, epsilon) {
  MLEs <- bind_cols(
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
    left_join(MLEs, by=c("cluster", "dim")) %>% 
    mutate(nll = -dgamma(value, shape=shape*(1-epsilon)/epsilon, rate=rate, log=TRUE)) %>% 
    filter(!is.infinite(nll)) %>%
    summarise(nll = sum(nll)) %>%
    pull
}
