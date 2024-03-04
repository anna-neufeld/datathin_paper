library(datathin)
setwd("~/datathin_paper/Figure2")

get_stepwise_split_cis <- function(train, test) {
  mod.0 <- lm(y~1, data=train)
  mod.full <- lm(y~., data=train)
  
  forward <- step(mod.0, direction='forward', scope=formula(mod.full), trace=0)
  
  mod.test <- lm(formula(forward), data=test)
  
  p <- NCOL(train)-1
  results <- matrix(NA, nrow=p, ncol=2)
  rownames(results) <- paste0("X",1:p, sep="")
  
  confints <- confint(mod.test)[-1,]
  results[rownames(confints),] <- confints
  
  return(results)
}

#### Run the sim
nTrials = 10000
n <- 200
p <- 20
p.imp <- 5

full.results.ss <- array(NA, dim=c(nTrials*3, p, 2)) 
full.results.datathin <-  array(NA, dim=c(nTrials*3, p, 2)) 



betas <- rep(rep(seq(0.01, 1, length.out=20), nTrials/20), each=3)
epsilons <- rep(c(0.2, 0.5, 0.8), nTrials)
sigma.sq <- 1

counter <- 1
for (t in 1:nTrials) {
  if (t%%100==1) {print(t)}
  set.seed(t)
  X <- matrix(rnorm(n*p), nrow=n, ncol=p)
  this.beta <- c(rep(betas[counter], p.imp), rep(0, p-p.imp))
  y <- X%*%this.beta+rnorm(n, mean=0, sd=sqrt(sigma.sq))
  dat <- data.frame(X,y)
  
  for (epsilon in c(0.2,0.5, 0.8)) {
    train.indices <- sample(1:n, size=epsilon*n)
    train.ss <- dat[train.indices, ]
    test.ss <- dat[-train.indices, ]
    full.results.ss[counter,,] <- get_stepwise_split_cis(train.ss,test.ss)
    
    thin.y <- datathin(y, family="normal", epsilon=c(epsilon, 1-epsilon), arg=sigma.sq)
    train.dt <- test.dt <- dat
    train.dt$y <- thin.y[,,1]
    test.dt$y <- thin.y[,,2]
    full.results.datathin[counter,,] <- get_stepwise_split_cis(train.dt,test.dt)
    counter <- counter + 1 
  }
}

power.check <- function(ints, truths) {
  return(ints[,1] > 0 | ints[,2] < 0)
}


### Check coverage
res.ss <- data.frame(t(sapply(1:(nTrials*3), function(u) power.check(full.results.ss[u,,], c(rep(betas[u], p.imp), rep(0, p-p.imp))))))
colnames(res.ss) <- paste0("X", 1:p, sep="") 
res.ss$beta <- betas
res.ss$epsilon <- epsilons
res.ss$method <- "Sample splitting"

res.dt <- data.frame(t(sapply(1:(nTrials*3), function(u) power.check(full.results.datathin[u,,], c(rep(betas[u], p.imp), rep(0, p-p.imp))))))
colnames(res.dt) <- paste0("X", 1:p, sep="") 
res.dt$beta <- betas
res.dt$epsilon <- epsilons
res.dt$method <- "Data thinning"

res <- rbind(res.ss, res.dt)

res.RANDOM <- rbind(res.ss, res.dt)

save(res.RANDOM, file="randomX_res-nov6.RData")
