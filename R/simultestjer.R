source('R/simul_1.R')
source('R/bcf_linear.R')

library(Rcpp)
sourceCpp("src/bcf_clean_overpar_linear.cpp")

library(ggplot2)
library(latex2exp)
library(rpart)
library(rpart.plot)
library(partykit)

n_simul <- 50
heter <- c(TRUE, FALSE)
linear <-c(TRUE, FALSE)
n  <- c(250, 500)
n_burn <- 200
n_sim <- 1000

data <- generate_data(500, is_te_hetero = TRUE, is_mu_nonlinear = TRUE, seed = 1995, RCT =T)
weights <- rep(1,500)
bcf_out <- bcf_linear(y                = data$y,
                      z                = data$z,
                      x_control        = as.matrix(data[, c(1:7)]),
                      x_moderate       = as.matrix(data[, c(1:7)]),
                      pihat            = data$pi_x,
                      nburn            = n_burn,
                      nsim             = n_sim,
                      w                = weights,
                      n_chains         = 2,
                      random_seed      = i,
                      update_interval  = 2000, 
                      no_output        = TRUE,
                      use_bscale = FALSE,
                      use_muscale = FALSE,
                      do_parallel = TRUE,
                      intTreat = TRUE)

load("C:/Users/P094412/Documents/bcf-linear/simulations output/bcf_out_het_TRUE_lin_FALSE_n_500_sim_50.RData")
library(coda)
chain_1 <- as.mcmc(bcf_out$beta)
effectiveSize(chain_1)
summary(bcf_out$beta_int)

hist(bcf_out$beta[,5])