
#auxiliary function to retrieve cred int from posterior
postsum <- function(posters,conf=0.05) {
  qs <- quantile(posters,p=c(conf/2,1-conf/2))
  mn <- mean(posters)
  return(c(qs[1],mn,qs[2]))
}

#intervals from posterior
intfrompost1 <- function(post,conf=0.05){
  shtotpost <- t(apply(post,1,postsum,conf=conf))
  ress <- list(tot=shtotpost)
  return(ress)
}

conftol <- function(margconf = 0.05, post, tol=0.05){
  #margconf <- 0.05; tol <- 0.1; post <- shapposts1
  int <- intfrompost1(post,conf=margconf)[[1]]
  difl <- ((post - int[,1]) >= 0) * ((post - int[,3]) <= 0)
  cs <- colMeans(difl)
  confret <- mean(cs >= 1-tol)
  return(confret)
} 

findconfglob <- function(post,tol=0.05){
  grid <- seq(0.0005, 0.05, by=0.0005)
  conftols <- sapply(grid, conftol, post=post, tol=tol)
  wm <- which(conftols - 0.95 < 0)[1] - 1
  confmarg <- grid[wm]
  return(confmarg)
}