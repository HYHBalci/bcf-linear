
all_beta <- c()
chains <- list()
for(iChain in 1:1) {
  # Create a dynamic filename based on iChain
  filename <- sprintf("chain_data_%d.RData", iChain)
  
  # Save the object to an RDS file
  load(filename)
  if(iChain==1){
    chain_1 <- chain_data
  } else if(iChain == 2) {
    chain_2 <- chain_data
  } else if(iChain == 3) {
    chain_3 <- chain_data
  } else {
    chain_4 <- chain_data
  }
  all_beta <- cbind(all_beta, t(chain_data$Beta))
}
library(coda)
chains <- mcmc.list(
  mcmc(as.matrix(chain_1$Beta[,3])),
  mcmc(as.matrix(chain_4$Beta[,3]))
)

gelman.diag(chains)

# all_beta <- t(all_beta)
chain_out <- vector("list", 4)

for (iChain in 1:4) {
  
  filename <- sprintf("chain_data_%d.RData", iChain)
  load(filename)
  # Run the BCF function and store the result
  chain_out[[iChain]] <- chain_data
  


}


all_sigma = c()
all_mu_scale = c()
all_tau_scale = c()
all_alpha = c()
all_b0 = c()
all_b1 = c()

all_yhat = c()
all_mu   = c()
all_tau  = c()
all_beta = c()
chain_list=list()

n_iter = length(chain_out[[1]]$sigma)

for (iChain in 1:4){
  
  sigma            <- chain_out[[iChain]]$sigma
  mu_scale         <- chain_out[[iChain]]$mu_scale
  tau_scale        <- chain_out[[iChain]]$tau_scale
  
  b0               <- chain_out[[iChain]]$b0
  b1               <- chain_out[[iChain]]$b1
  
  beta <- chain_out[[iChain]]$Beta
  alpha <- chain_out[[iChain]]$alpha
  
  yhat             <- chain_out[[iChain]]$yhat
  tau              <- chain_out[[iChain]]$tau
  mu               <- chain_out[[iChain]]$mu
  has_file_output  <- chain_out[[iChain]]$has_file_output
  
  # -----------------------------    
  # Support Old Output
  # -----------------------------
  all_sigma       = c(all_sigma,     sigma)
  all_mu_scale    = c(all_mu_scale,  mu_scale)
  all_tau_scale   = c(all_tau_scale, tau_scale)
  all_b0 = c(all_b0, b0)
  all_b1 = c(all_b1, b1)
  
  all_yhat = rbind(all_yhat, yhat)
  all_mu   = rbind(all_mu,   mu)
  all_tau  = rbind(all_tau,  tau)
  
  all_beta = rbind(all_beta, beta)
  all_alpha = rbind(all_alpha, alpha)
  # -----------------------------    
  # Make the MCMC Object
  # -----------------------------
  w <- matrix(1, ncol = 1, nrow = 1000)
  scalar_df <- data.frame("sigma"     = sigma,
                          "b0"  = b0, 
                          "b1"  = b1)
  
  # y_df <- as.data.frame(chain$yhat)
  # colnames(y_df) <- paste0('y',1:ncol(y_df))
  # 
  # mu_df <- as.data.frame(chain$mu)
  # colnames(mu_df) <- paste0('mu',1:ncol(mu_df))
  # 
  # tau_df <- as.data.frame(chain$tau)
  # colnames(tau_df) <- paste0('tau',1:ncol(tau_df))
  
  chain_list[[iChain]] <- coda::as.mcmc(scalar_df)
  # -----------------------------    
  # Sanity Check Constants Across Chains
  # -----------------------------
  # if(chain_out[[iChain]]$sdy              != chain_out[[1]]$sdy)              stop("sdy not consistent between chains for no reason")
#   if(chain_out[[iChain]]$con_sd           != chain_out[[1]]$con_sd)           stop("con_sd not consistent between chains for no reason")
#   if(chain_out[[iChain]]$mod_sd           != chain_out[[1]]$mod_sd)           stop("mod_sd not consistent between chains for no reason")
#   if(chain_out[[iChain]]$muy              != chain_out[[1]]$muy)              stop("muy not consistent between chains for no reason")
#   if(chain_out[[iChain]]$include_pi       != chain_out[[1]]$include_pi)       stop("include_pi not consistent between chains for no reason")
#   if(any(chain_out[[iChain]]$perm         != chain_out[[1]]$perm))            stop("perm not consistent between chains for no reason")
#   if(chain_out[[iChain]]$has_file_output  != chain_out[[1]]$has_file_output)  stop("has_file_output not consistent between chains for no reason")
  }

fitObj <- list(sigma = all_sigma,
               yhat = all_yhat,
               sdy = chain_out[[1]]$sdy,
               muy = chain_out[[1]]$muy,
               mu  = all_mu,
               tau = all_tau,
               beta = all_beta,
               alpha = all_alpha,
               mu_scale = all_mu_scale,
               tau_scale = all_tau_scale,
               b0 = all_b0,
               b1 = all_b1,
               include_pi = chain_out[[1]]$include_pi,
               random_seed = chain_out[[1]]$random_seed,
               coda_chains = coda::as.mcmc.list(chain_list),
               raw_chains = chain_out, 
               has_file_output = has_file_output)

attr(fitObj, "class") <- "bcf"

.cleanup_after_par(do_type_config)

windows()  # For Windows
hist(chain_data$m_post)


plot(sd(y)*chain_data$Beta[,3], type = "l", col = "blue", main = "Trace Plot of MCMC Draws",
     xlab = "Iteration", ylab = "Value")
abline(h = mean(chain_data$Beta[,3]), col = "red", lwd = 2, lty = 2) 

fitObj$beta_adj <- fitObj$beta*sd(y)

library(coda)
geweke.diag(chain_data$Beta[,2])
raftery.diag(chain_data$Beta[,2])