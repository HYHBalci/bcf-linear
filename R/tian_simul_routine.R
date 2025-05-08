source('R/simul_1.R')
source('R/test_linked_shrinkage.R')

n_simul <- 50
heter <- c(TRUE, FALSE)
linear <- c(TRUE,FALSE)
n <- c(250,500)

for (het in heter) {
  for (lin in linear) {
    for (n_obser in n) {
      for (i in 1:n_simul) {
        set.seed(i)
      
        data <- generate_data_2(n_obser, is_te_hetero = het, is_mu_nonlinear = lin, seed = i, RCT = TRUE, z_diff = T, tian = T)
        
        X <- as.matrix(data[, c("x1", "x2", "x3", "x4", "x5_1", "x5_2")])
        y <- 2 * data$y * data$z 
        z <- data$z
        
        # --- Run Tian Model ---
        posterior <- sample_linear_part(y, rep(1,n_obser), X, intTreat = TRUE, iter = 8000, burnin = 1000)
        
        filename <- sprintf("tian_heter_%s_linear_%s_n_%d_sim_%d.RData", 
                            ifelse(het, "T", "F"), 
                            ifelse(lin, "T", "F"), 
                            n_obser, 
                            i)
        save(posterior, file = filename)
        print(filename)
      }
    }
  }
}