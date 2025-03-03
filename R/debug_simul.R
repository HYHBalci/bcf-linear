source('R/simul_1.R')
source('R/bcf_linear.R')

library(Rcpp)
sourceCpp("src/bcf_clean_overpar_linear.cpp")

library(ggplot2)
library(latex2exp)
library(rpart)
library(rpart.plot)
library(partykit)
n_burn <- 2000
n_sim <- 5000
data <- generate_data(250, is_te_hetero = T, is_mu_nonlinear = T, seed = 1848, test = T, RCT = T)
weights <- rep(1,250)
bcf_out <- bcf_linear(y                = data$y,
                      z                = data$z,
                      x_control        = as.matrix(data[, c(1:7)]),
                      x_moderate       = as.matrix(data[, c(1:7)]),
                      pihat            = data$pi_x,
                      nburn            = n_burn,
                      nsim             = n_sim,
                      nthin = 1,
                      w                = weights,
                      n_chains         = 2,
                      random_seed      = i,
                      update_interval  = 100, 
                      no_output        = TRUE,
                      use_bscale = FALSE,
                      use_muscale = FALSE,
                      do_parallel = TRUE,
                      intTreat = TRUE,
                      hamiltonian = F)

load("simulations output/bcf_out_het_TRUE_lin_TRUE_n_250_sim_1.RData")

