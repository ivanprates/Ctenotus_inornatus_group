#############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### April 2023.
### The goals of this script are:
### To help define gamma-distributed priors parameters in G-PhoCS analyses.

# Given a mutation rate of:
u = 7.6*10^-9

## 1. For tau:
## Tree height (in number of generations): 
gen <- 5000
gen <- 500000
r = u*gen 
r ## From 0.00004 to 0.004.

## 2. For theta:
## Effective population size:
N <- 100000
N <- 1000000
theta <- 4*N*u
theta ## From 0.0003 to 0.03.

## 3. For migration rate:
## M = proportion of ind. in pop. 1 that arose from migration from pop. 2 per generation.

## m for different M:
M <- 1/200000 ## Under N = 100,000, this leads to 2NM = 1.
M <- 1/2000000 ## Under N = 1,000,000, this leads to 2NM = 1.
m = M/u
m ## from 66 to 658.
2*N*M

## To achieve these ranges, implement gamma parameters as:
shape <- 1 ; rate <- 100 ; xlmax <- 0.05 ## For tau and theta.
shape <- 0.001 ; rate <- 0.00001 ; xlmax <- 1000 ## For migration.

## Plot:
xmax <- qgamma(p = 0.999, shape = shape, rate = rate)
x <- seq(from = 0, to = xmax, by = xmax/1000)
dens <- dgamma(x, shape = shape, rate = rate)
plot(x, dens, type = 'l', lwd = 2, xlab = "parameter", xlim = c(0, xlmax))

## End of script.
