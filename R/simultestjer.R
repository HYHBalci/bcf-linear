source('R/simul_1.R')
source('R/bcf_linear.R')

library(Rcpp)
sourceCpp("src/bcf_clean_overpar_linear.cpp")

n_simul <- 50
heter <- c(TRUE, FALSE)
linear <-c(TRUE, FALSE)
n  <- c(250, 500)
n_burn <- 500
n_sim <- 3000
counter = 1


library(ggplot2)
library(latex2exp)
library(rpart)
library(rpart.plot)
library(partykit)

n_simul <- 50
heter <- c(TRUE, FALSE)
linear <-c(TRUE, FALSE)
n  <- c(250, 500)
n_burn <- 1000
n_sim <- 2000
counter = 1
n_obser = 500
data <- generate_data(n_obser, is_te_hetero = TRUE, is_mu_nonlinear = F, seed = 18)
weights <- rep(1,n_obser)
bcf_out <- bcf_linear(y                = data$y,
                      z                = data$z,
                      x_control        = as.matrix(data[, c(1:7)]),
                      x_moderate       = as.matrix(data[, c(1:7)]),
                      pihat            = data$pi_x,
                      nburn            = n_burn,
                      nsim             = n_sim,
                      nthin = 1,
                      ntree_control = 400,
                      w                = weights,
                      n_chains         = 2,
                      random_seed      = i,
                      update_interval  = 2000, 
                      no_output        = TRUE,
                      use_bscale = FALSE,
                      use_muscale = FALSE,
                      do_parallel = TRUE,
                      intTreat = TRUE, hamiltonian = F, sparse = T)

load('2x10000-sparse.RData')
chain_1 <- as.mcmc(bcf_out$raw_chains[[1]]$beta_int)
chain_2 <- as.mcmc(bcf_out$raw_chains[[2]]$Beta)
mcmc_list <- mcmc.list(chain_1, chain_2)
gelman.diag(mcmc_list)
chain_sigma <- as.mcmc(bcf_out$raw_chains[[1]]$sigma)
eff_300 <- effectiveSize(chain_1)
traceplot(chain_1)
hist(bcf_out$alpha*sd(data$y) + mean(data$y))
mean(bcf_out$alpha*sd(data$y) + mean(data$y))
hist(bcf_out$alpha*bcf_out$sdy + bcf_out$muy, breaks = 30)
save(bcf_out, file = '2x10000-sparse-false.RData')
bcf_out$muy
bcf_out$sdy


hist(bcf_out$beta_int[,2], breaks = 200)