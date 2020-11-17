# Problem 1
## 1a: Histograms

d <- read.csv("http://www.stat.yale.edu/~jtc5/238/data/cost-of-the-muse.csv")
novels <- d[d$Type1 == 1, 3]
poems <- d[d$Type1 == 2, 3]
nonfiction <- d[d$Type1 == 3, 3]

par(mfrow=c(3, 1))
rg <- range(d[, 3])
hist(novels, 100,  col="red", xlim=rg)
hist(poems, 100,  col="red", xlim=rg)
hist(nonfiction, 100,  col="red", xlim=rg)

## 1b: MCMC and Quantiles
lik <- function(th) {
  mu1 <- th[1]; sig1 <- th[2]; mu2 <- th[3]; sig2 <- th[4]; mu3 <- th[5]; sig3 <- th[6]
  return(
    prod(dnorm(x=novels, mean=mu1, sd=sig1)) *
      prod(dnorm(x=poems, mean=mu2, sd=sig2)) *
      prod(dnorm(x=nonfiction, mean=mu3, sd=sig3))
  )
}

prior <- function(th) {
  mu1 <- th[1]; sig1 <- th[2]; mu2 <- th[3]; sig2 <- th[4]; mu3 <- th[5]; sig3 <- th[6]
  return(
    dunif(mu1, min = 0, max = 100) *
      dunif(mu2, min = 0, max = 100) *  
      dunif(mu3, min = 0, max = 100) * 
      dunif(sig1, min = 0, max = 100) * 
      dunif(sig2, min = 0, max = 100) *
      dunif(sig3, min = 0, max = 100)
  )
}

post <- function(th) {
  if((th[2] < 0) | (th[4] < 0) | (th[6] < 0)) {
    return(0)
  }
  
  return(prior(th) * lik(th))
}

nit <- 100000
results <- matrix(0, nrow = nit, ncol = 6)
th <- c(70, 10, 70, 10, 70, 10)   # Initial value
results[1, ] <- th
for(i in 2:nit) {
  cand <- th + rnorm(6)
  ratio <- post(cand)/post(th)
  lik(th)
  u <- runif(n = 1, min = 0, max = 1)
  if(u < ratio){
    th <- cand
  }
  results[i,] <- th
}

r <- data.frame(results)
names(r) <- c("mu1", "sig1", "mu2", "sig2", "mu3", "sig3")

quantile(r[,1], prob=c(0.025, 0.5, 0.975))
quantile(r[,2], prob=c(0.025, 0.5, 0.975))
quantile(r[,3], prob=c(0.025, 0.5, 0.975))
quantile(r[,4], prob=c(0.025, 0.5, 0.975))
quantile(r[,5], prob=c(0.025, 0.5, 0.975))
quantile(r[,6], prob=c(0.025, 0.5, 0.975))

## 1c: Deltas
deltas1 <- r[,1] - r[,3]
deltas2 <- r[,5] - r[,3]
hist(deltas1, col=5, prob=T)
hist(deltas2, col=5, prob=T)
quantile(deltas1, prob=c(0.025, 0.975))
quantile(deltas2, prob=c(0.025, 0.975))

## 1d: Hierarchies
# Probability of mu2 < mu1 < mu3
sum(r[,3] < r[,1] & r[,1] < r[,5])/nit

# Other orderings:
sum(r[,1] < r[,3] & r[,3] < r[,5])/nit
sum(r[,1] < r[,5] & r[,5] < r[,3])/nit
sum(r[,3] < r[,5] & r[,5] < r[,1])/nit
sum(r[,5] < r[,1] & r[,1] < r[,3])/nit
sum(r[,5] < r[,3] & r[,3] < r[,1])/nit

# Problem 2: JAGS Implementation
library(rjags)
mymodel <- "
  model{
    for(i in 1:67) {
      novels[i] ~ dnorm(mu1, tau1)
    }
    for(i in 1:32) {
      poems[i] ~ dnorm(mu2, tau2)
    }
    for(i in 1:24) {
      nonfiction[i] ~ dnorm(mu3, tau3)
    }
    
    mu1 ~ dunif(0, 100)
    mu2 ~ dunif(0, 100)
    mu3 ~ dunif(0, 100)
    sig1 ~ dunif(0, 100)
    sig2 ~ dunif(0, 100)
    sig3 ~ dunif(0, 100)
    tau1 <- 1/(sig1^2)
    tau2 <- 1/(sig2^2)
    tau3 <- 1/(sig3^2)
  }
"

jm <- jags.model(file=textConnection(mymodel),
                 data=list(novels=novels, poems=poems, nonfiction=nonfiction),
                 inits=list(mu1=50, sig1=10, mu2=50, sig2=10, mu3=50, sig3=10))
cs <- coda.samples(jm, c("mu1", "sig1", "mu2", "sig2", "mu3", "sig3"), nit)
s <- as.data.frame(cs[[1]])

quantile(s[,1], prob=c(0.025, 0.5, 0.975))
quantile(s[,2], prob=c(0.025, 0.5, 0.975))
quantile(s[,3], prob=c(0.025, 0.5, 0.975))
quantile(s[,4], prob=c(0.025, 0.5, 0.975))
quantile(s[,5], prob=c(0.025, 0.5, 0.975))
quantile(s[,6], prob=c(0.025, 0.5, 0.975))

## 2b
deltas1 <- s[,1] - s[,2]
deltas2 <- s[,3] - s[,2]
hist(deltas1, col=5, prob=T)
hist(deltas2, col=5, prob=T)
quantile(deltas1, prob=c(0.025, 0.975))
quantile(deltas2, prob=c(0.025, 0.975))

## 2c
# Probability of mu2 < mu1 < mu3
sum(s[,2] < s[,1] & s[,1] < s[,3])/nit

# Other orderings:
sum(s[,1] < s[,2] & s[,2] < s[,3])/nit
sum(s[,1] < s[,3] & s[,3] < s[,2])/nit
sum(s[,2] < s[,3] & s[,3] < s[,1])/nit
sum(s[,3] < s[,1] & s[,1] < s[,2])/nit
sum(s[,3] < s[,2] & s[,2] < s[,1])/nit


# Problem 3: Helsinki Heart Study
## 3a: JAGS Implementation
treatments <- 56
controls <- 84
mymodel <- "
  model {
    treatments ~ dbin(th_treatment, 1000)
    controls ~ dbin(th_control, 1000)
    th_treatment ~ dunif(0, 1)
    th_control ~ dunif(0, 1)
  }
"

jm <- jags.model(textConnection(mymodel), data=list(treatments=treatments, controls=controls))
cs <- coda.samples(jm, c("th_treatment", "th_control"), 10000)
s <- as.data.frame(x=cs[[1]])
hist(s$th_treatment, xlim=c(0, 0.5))
hist(s$th_control, xlim=c(0, 0.5))

## 3b: Quantiles
quantile(s[,1], prob=c(0.025, 0.975))
quantile(s[,2], prob=c(0.025, 0.975))

## 3c: Reduction
reduction = ((s[,1] - s[,2]) * 100)/s[,1]
quantile(reduction, prob=c(0.025, 0.975))

# Problem 4: Implementation using Conjugate Prior (Beta)
qbeta(p=c(.025, .975), shape1=85, shape2=917)
qbeta(p=c(.025, .975), shape1=57, shape2=945)
sample1 <- rbeta(10000, shape1=85, shape2=917)
sample2 <- rbeta(10000, shape1=57, shape2=945)
reduction = ((sample1 - sample2) * 100)/sample1
quantile(reduction, prob=c(0.025, 0.975))