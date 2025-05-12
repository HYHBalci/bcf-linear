# dummy versus deviance coding
n <- 100
abc <- sample(c("a", "b", "c"), n, replace = TRUE)
y <- 1 + 1 * I(abc=="a") - 2 * I(abc=="b") + 1 * I(abc=="c") + rnorm(n, 0, .1) 
abc <- factor(abc)
contrasts(abc) # default dummy coding against lowest reference category
summary(lm(y ~ 1 + abc))
# true coefs are 2 (int + beta_a), -3 (difference beta_b and beta_a) and 0 (difference beta_c and beta_a)
contr.sum(3) # deviance coding
contrasts(abc) <- contr.sum(3) 
mod <- lm(y ~ 1 + abc)
summary(mod)
# get all coefficients for a,b,c
contr.sum(3) %*% coef(mod)[-1]
# overall average is indeed 1 (now the intercept)
# a vs overall equals +1
# b vs overall equals -2
# c vs overall is +1
