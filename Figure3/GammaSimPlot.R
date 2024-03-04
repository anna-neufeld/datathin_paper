library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(1234)

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

# Generate Gamma cluster data
## Number of data points per cluster
n <- 100
## Number of dimensions
d <- 2
## Number of clusters
ncluster <- 4

## Gamma parameters for the "small d, small K" setting
pars <- data.frame(lambda = c(20,20,20,20,20,20,20,20), 
                   theta = c(0.5,5,5,0.5,10,10,0.5,0.5))
lambda <- matrix(rep(pars$lambda,n), nrow=ncluster*n, ncol=d, byrow=TRUE)
theta <- matrix(rep(pars$theta,n), nrow=ncluster*n, ncol=d, byrow=TRUE)

## Generate the data
X <- matrix(0, nrow=ncluster*n, ncol=d)
X[] <- rgamma(n*d*ncluster, shape=lambda, rate=theta)

# Data thin the Gamma data with epsilon=0.5
Xdt <- gammasplit(X, epsilon=0.5, shape=lambda)
Xall <- bind_rows(
  as_tibble(X) %>% mutate(dataset="X"),
  as_tibble(purrr::map_dfr(Xdt, ~as.data.frame(.x), .id="dataset"))
) %>% 
  mutate(dataset = case_when(
    dataset == "Xte" ~ "X^{(2)}",
    dataset == "Xtr" ~ "X^{(1)}",
    TRUE ~ "X"
  ))


# Create and display the plot
p <- Xall %>% 
  mutate(g = as.factor(rep(c(1:ncluster), n*3))) %>%
  ggplot(aes(x=V1, y=V2, colour=g)) +
  geom_point() +
  facet_wrap(~dataset, labeller = label_parsed) +
  xlab("") + ylab("") +
  theme_bw() + 
  theme(legend.position = "none")

p+theme(axis.text=element_text(size=16), strip.text=element_text(size=16))
setwd("~/Dropbox/Generalized Count Splitting/JMLR-resubmit-sep-2023/Figures/")
ggsave("2dgammasim-1.png", width=10, height=4)
