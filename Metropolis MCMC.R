# Problem 2
## 2a: Example of POOR Scale Size
f <- function(x) {
  return(1/(3*sqrt(2*pi)) * exp(-0.5*((x - 3 )/ 0.001)^2))
}

nit <- 10000
path <- rep(0, nit)
state <- 3  # Initial state
path[1] <- state
scale <- 1
for (i in 2:nit) {
  candidate <- runif(n = 1, min = state - scale, max = state + scale)
  ratio <- f(candidate) / f(state)
  u <- runif(n = 1, min = 0, max = 1)
  if (u < ratio) {
    state <- candidate
  }
  path[i] <- state
}

plot(path[-(1:200)])

library(MASS)
truehist(path[-(1:200)])
con <- integrate(f, 2.9, 3.1)$value     
fnormalized <- function(x){f(x)/con}
curve(fnormalized(x), add=TRUE, col="red", lwd=2)

## 2b: Example of GOOD Scale size
path <- rep(0, nit)

state <- 3  # Initial state
path[1] <- state

scale <- 0.001
for (i in 2:nit) {
  candidate <- runif(n = 1, min = state - scale, max = state + scale)
  ratio <- f(candidate) / f(state)
  u <- runif(n = 1, min = 0, max = 1)
  if (u < ratio) {
    state <- candidate
  }
  path[i] <- state
}

plot(path[-(1:200)])

library(MASS)
truehist(path[-(1:200)])
curve(fnormalized(x), add=TRUE, col="red", lwd=2)

## 2c: Another example of POOR scale size
f <- function(x) {
  return(1/sqrt(2*pi) * exp(-(x^2)/2))
}

path <- rep(0, nit)

state <- 3  # Initial state
path[1] <- state

scale <- 0.001
for (i in 2:nit) {
  candidate <- runif(n = 1, min = state - scale, max = state + scale)
  ratio <- f(candidate) / f(state)
  u <- runif(n = 1, min = 0, max = 1)
  if (u < ratio) {
    state <- candidate
  }
  path[i] <- state
}

plot(path[-(1:200)])
con <- integrate(f, -5, 5)$value     
fnormalized <- function(x){f(x)/con}

library(MASS)
truehist(path[-(1:200)])
curve(fnormalized(x), add=TRUE, col="red", lwd=2)

# Problem 4
## 4b: Independence sampler MCMC
f <- function(x) {
  return((1 + (x^2)/3)^(-0.5*4))
}

nit <- 100000
path <- rep(NA, nit)
state <- 0 # Initial state
path[1] <- state
for (i in 2:nit) {
  angle <- runif(n = 1, min = -pi/2, max = pi/2)
  candidate <- tan(angle)
  ratio <- (f(candidate) * dcauchy(state)) / (f(state) * dcauchy(candidate))
  u <- runif(n = 1, min = 0, max = 1)
  if (u < ratio) {
    state <- candidate # Accept the candidate
  }
  path[i] <- state
}

plot(path[-(1:200)])

## 4c: Quantiles
quantile(abs(path), probs = c(0.9, 0.95, 0.99))
# Factoring in the absolute values, for comparison:
qt(p = c(0.95, 0.975, 0.995), df = 3)

# Problem 5
## 5b: Hardcore Lattice model
nit <- 100000
path <- array(0, dim=c(8, 8, nit))

meets_constraint <- function(lattice) {
  # Check for vertical pairs
  for (i in 1:8) {
    for (j in 1:7) {
      if (lattice[i, j] + lattice[i, j + 1] == 2) {
        return(FALSE)
      }
    }
  }
  
  # Check for horizontal pairs
  for (j in 1:8) {
    for (i in 1:7) {
      if (lattice[i, j] + lattice[i + 1, j] == 2) {
        return(FALSE)
      }
    }
  }
  
  return(TRUE);
}

for (k in 2:nit) {
  x <- sample(1:8, size=1)
  y <- sample(1:8, size=1)
  this_lattice <- path[,,k - 1]
  this_lattice[x, y] <- 1 - this_lattice[x, y]
  if (meets_constraint(this_lattice)) {
    path[,,k] <- this_lattice
  }
  else {
    path[,,k] <- path[,,k - 1]
  }
}

occupancies <- rep(0, nit)
for (i in 1:nit) {
  occupancies[i] <- sum(path[,,i])
}

plot(table(occupancies))

## 5c: Quantiles
quantile(occupancies, c(0.1, 0.5, 0.9))

## 5d: Graphs
library(grid)

plot.hc <- function(mat){
  m <- dim(mat)[1]
  n <- dim(mat)[2]
  
  grid.newpage()
  rad = .02
  xx = function(j) {j/(n+1)}
  yy = function(i) {(m+1-i)/(m+1)}
  
  for(i in 1:m){
    grid.lines(c(xx(1),xx(n)),c(yy(i),yy(i)))
  }
  for(j in 1:n){
    grid.lines(c(xx(j),xx(j)),c(yy(1),yy(m)))
  }
  for(j in 1:n){
    for(i in 1:m){
      mycolor <- c('black','green')[1+mat[i,j]]
      grid.circle(xx(j), yy(i), rad, gp=gpar(fill=mycolor))
    }
  }
}

plot.hc(path[,,2500])
plot.hc(path[,,5000])
plot.hc(path[,,7500])
plot.hc(path[,,10000])