source('R/bcf_linear.R')
library(Rcpp)
sourceCpp("src/bcf_clean_overpar_linear.cpp")
sourceCpp("src/bcf_linear_2.cpp")

library(bcf)
library(ggplot2)
library(latex2exp)
library(rpart)
library(rpart.plot)
library(partykit)

set.seed(1)

p <- 3 # two control variables and one effect moderator
n <- 1000 # number of observations
n_burn <- 2000
n_sim <- 5000

x <- matrix(rnorm(n*(p)), nrow=n)
# x <- cbind(x, x[,2] + rnorm(n))
weights <- rep(1,n)


# create targeted selection, whereby a practice's likelihood of joining the intervention (pi) 
# is related to their expected outcome (mu)
mu <- -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2])) - 0.1
# generate treatment variable
pi <- rep(0.5,n)
z <- rbinom(n,1,pi)

# tau is the true treatment effect. It varies across practices as a function of
# X3 and X2 <- 2*x[,3] 
tau <- 0.1*x[,2] + 2*x[,3] + 0.5*x[,1]

# generate the expected response using mu, tau and z
y_noiseless <- mu + tau*z

# set the noise level relative to the expected mean function of Y
sigma <- 1
set.seed(1848)
# draw the response variable with additive error
y <- y_noiseless + sigma*rnorm(n)
bcf_out <- bcf_linear(y                = y,
                      z                = z,
                      x_control        = x,
                      x_moderate       = x,
                      pihat            = pi,
                      nburn            = n_burn,
                      nsim             = n_sim,
                      w                = weights,
                      n_chains         = 4,
                      random_seed      = 1,
                      update_interval  = 1, 
                      no_output        = FALSE,
                      use_bscale = FALSE,
                      use_muscale = FALSE,
                      do_parallel = TRUE,
                      tian = F,
                      test = F)