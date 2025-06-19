library(coda)

# Load and compute ESS
ess_results <- data.frame(Step_Index = numeric(), ESS = numeric())

for (i in 1:length(step_outs)) {
  filename <- sprintf("BCF_fit_step_sim_%d.Rdata", i)
  load(filename)
  
  # Assume nbcf_fit$beta_samples has dimension (chains, iterations, parameters)
  beta_samples <- nbcf_fit$Beta
  
  # Combine chains
  mcmc_list <- lapply(1:dim(beta_samples)[1], function(chain) {
    mcmc(beta_samples[chain,,], start=1)
  })
  
  mcmc_combined <- mcmc.list(mcmc_list)
  ess_vals <- effectiveSize(mcmc_combined)
  
  ess_results <- rbind(ess_results, data.frame(Step_Index = i, ESS = mean(ess_vals)))
}

print(ess_results)
