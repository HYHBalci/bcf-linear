source('R/bcf_linear.R')
library(Rcpp)
sourceCpp("src/bcf_clean_overpar_linear.cpp")
library(bcf)
library(ggplot2)
library(latex2exp)
library(rpart)
library(rpart.plot)
library(partykit)

set.seed(1)

p <- 3 # two control variables and one effect moderator
n <- 500 # number of observations
n_burn <- 300
n_sim <- 1000

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
# X3 and X2
tau <- 2*x[,3] + 0.1*x[,2] + x[,2]*x[,3]

# generate the expected response using mu, tau and z
y_noiseless <- mu + tau*z

# set the noise level relative to the expected mean function of Y
sigma <- 1

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
               n_chains         = 1,
               random_seed      = 1,
               update_interval  = 1, 
               no_output        = FALSE,
               use_bscale = FALSE,
               use_muscale = FALSE,
               do_parallel = TRUE,
               intTreat = TRUE)

bcf_out$beta

summary(bcf_out)


tau_ests <- data.frame(Mean  = colMeans(bcf_out$b_p),
                       Low95 = apply(bcf_out$tau, 2, function(x) quantile(x, 0.025)),
                       Up95  = apply(bcf_out$tau, 2, function(x) quantile(x, 0.975)))

ggplot(NULL, aes(x = x[,2])) +
  geom_pointrange(aes(y = tau_ests$Mean, ymin = tau_ests$Low95, ymax = tau_ests$Up95), color = "forestgreen", alpha = 0.5) +
  geom_smooth(aes(y = tau), se = FALSE) +
  xlab(TeX("$x_3$")) +
  ylab(TeX("$\\hat{\\tau}$")) +
  xlim(-4, 6) +
  geom_segment(aes(x = 3, xend = 4, y = 0.2, yend = 0.2), color = "blue", alpha = 0.9) +
  geom_text(aes(x = 4.5, y = 0.2, label = "Truth"), color = "black") +
  geom_segment(aes(x = 3, xend = 4, y = 0.1, yend = 0.1), color = "forestgreen", alpha = 0.7) +
  geom_text(aes(x = 5.2, y = 0.1, label = "Estimate (95% CI)"), color = "black")
#> `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'
#> Warning: Removed 6 rows containing non-finite values (`stat_smooth()`).
#> Warning: Removed 6 rows containing missing values (`geom_pointrange()`).
#> 
#> 
#> 
#> 1) Tian Solution 
#> 2) Mu_post z_i correlation
#> 3) Tree structure, sparse prior
