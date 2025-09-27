#' Bayesian Adaptive Polynomial Chaos Expansions
#'  This code reproduces the simulation study used in the paper.
#'  This doesn't need to be run, as the data is stored in the
#'  data subdirectory. See duqling package documentation for
#'  details (github.com/knrumsey/duqling)

library(duqling) # devtools::install_github("knrumsey/duqling")
library(khaos)   # devtools::install_github("knrumsey/khaos")
library(BART)
library(laGP)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stargazer) # For tables, not super necesarry


#=================================
#  EMULATOR SPECIFIC FUNCTIONS
#=================================
fit1 <- function(X, y){
  adaptive_khaos(X, y, legacy=TRUE)
}
fit2 <- function(X, y){
  adaptive_khaos2(X, y)
}
fit3 <- function(X, y){
  sparse_khaos(X, y, max_basis=1e5)
}
fit4 <- function(X, y){
  BART::wbart(X, y)
}
fit5 <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-6)
  if(max(table(y)) > 10) y <- y + rnorm(y, 0, 1e-6)
  list(X=X, y=y)
}
pred5 <- function(obj, Xt){
  n <- nrow(Xt)
  nn <- min(100, max(30, floor(sqrt(n))))

  preds <- matrix(NA, nrow=1000, ncol=nrow(Xt))
  for(i in 1:nrow(Xt)){
    mod <- laGP(Xt[i,,drop=FALSE], start=10, end=nn, X=obj$X, Z=obj$y)
    m <- mod$mean
    v <- mod$df
    s2 <- mod$s2 * (v-2) / v
    preds[,i] <- m + rt(1000, v) * sqrt(s2)
  }
  return(preds)
}
pred <- function(obj, Xt){
  predict(obj, Xt)
}
ff <- list(fit1, fit2, fit3, fit4, fit5)
pp <- list(pred, pred, pred, pred, pred5)


#=================================
#  SIM STUDY WITH DUQLING
#=================================
if(rerun_simulations){
  res <- duqling::run_sim_study(ff, pp,
                                fnames=c("rabbits", "ishigami",
                                         "pollutant_uni", "banana",
                                         "friedman20"),
                                n_train=1000, NSR=c(0, 0.5), replications=10,
                                mc_cores=10,
                                method_names = c('KHAOS (ridge)',
                                                 'KHAOS (g-prior)',
                                                 'sparsePCE',
                                                 'BART',
                                                 'laGP'))

}
