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
counter = 1
for(het in heter){
  for(lin in linear){
    for(n_obser in n){
      for(i in 1:n_simul){
        data <- generate_data(n_obser, is_te_hetero = het, is_mu_nonlinear = lin, seed = i)
        weights <- rep(1,n_obser)
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
        # Define the file name for saving the result
        file_name <- paste0("bcf_out_het_", het, "_lin_", lin, "_n_", n_obser, "_sim_", i, ".RData")
        
        # Save the bcf object
        save(bcf_out, file = file_name)
        flush.console()
      }
    }
  }
}