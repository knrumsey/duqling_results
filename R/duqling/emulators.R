library(duqling)
library(laGP)

# MASTER FUNCTIONS
fit_lagp <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-6)
  if(max(table(y)) > 10) y <- y + rnorm(y, 0, 1e-6)
  list(X=X, y=y)
}
pred_lagp <- function(obj, Xt, nn, mm){
  n <- nrow(Xt)

  preds <- matrix(NA, nrow=1000, ncol=nrow(Xt))
  for(i in 1:nrow(Xt)){
    nug <- max(1e-7 * var(obj$y), 1e-9)
    mod <- laGP(Xt[i,,drop=FALSE], start=10, end=nn, X=obj$X, Z=obj$y, method=mm, g=nug)
    m <- mod$mean
    v <- mod$df
    s2 <- mod$s2 * (v-2) / v
    preds[,i] <- m + rt(1000, v) * sqrt(s2)
  }
  return(preds)
}


method_list <- c("alcray", "alc", "mspe", "nn")
nn_list <- c(25, 50, 100, 200)
method_names <- NULL
for (i in seq_along(method_list)) {
  for (j in seq_along(nn_list)) {
    mn_curr <- paste0(method_list[i], "_", nn_list[j])
    method_names <- c(method_names, mn_curr)

    # Use local to bind i and j
    assign(paste0("pred_", mn_curr),
           local({
             mm <- method_list[i]
             nn <- nn_list[j]
             function(obj, Xt) pred_lagp(obj, Xt, nn = nn, mm = mm)
           })
    )

    # Fit function is always the same
    assign(paste0("fit_", mn_curr), fit_lagp)
  }
}
model_names <- method_names
(base) [knrumsey@darwin-fe1 duqling_sims]$
  (base) [knrumsey@darwin-fe1 duqling_sims]$ ls
bnns                         check_failures_multi.R   emulators.R       out           rashomon.R              run_single_case.R            submit_single_cases.sh
bnns_0.1.2.tar.gz            check_lagp_missing.R     fname_metrics     out_data      rashomon_submit.sh      slurm_logs
ccs6_submit_single_cases.sh  custom_data_files.Rda    get_missing.R     out_lagp      run_single_case_data.R  submit_single_cases_data.sh
check_failures_multi_data.R  custom_duqling_data.Rda  lagp_emulators.R  rashomon.csv  run_single_case_lagp.R  submit_single_cases_lagp.sh
(base) [knrumsey@darwin-fe1 duqling_sims]$ cat emulators.R
# Save me some time during debugging:
library(duqling)
X <- lhs::randomLHS(100, 3)
y <- apply(X,1,duqling::ishigami)

# Keep track of model names and use get(paste0("fit_", model_name))

model_names <- NULL

# 1. BASS
library(BASS)
fit_bass <- function(X, y) bass(X, y, verbose=FALSE)
pred_bass <- function(obj, Xt) predict(obj, Xt)
model_names <- c(model_names, "bass")

# 2. BART
library(BART)
fit_bart <- function(X, y){
  if(sd(y) ==0) y <- y + rnorm(y,0,1e-7)
  wbart(X, y)
}
pred_bart <- function(obj, Xt) predict(obj, Xt)
model_names <- c(model_names, "bart")

# 3. BayesPPR
library(BayesPPR)
fit_bppr <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-6)
  bppr(X, y)
}
pred_bppr <- function(obj, Xt) predict(obj, Xt)
model_names <- c(model_names, "bppr")

# 4. Bayesian linear model
fit_blm <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-7)
  ytrain <- y
  Xtrain <- cbind(1, X)
  n <- nrow(Xtrain)
  p <- ncol(Xtrain)
  XtX_inv <- solve(t(Xtrain) %*% Xtrain)
  beta_hat <- XtX_inv %*% t(Xtrain) %*% ytrain
  residuals <- ytrain - Xtrain %*% beta_hat
  sigma2_hat <- as.numeric((t(residuals) %*% residuals) / (n - p))
  list(
    beta_hat = beta_hat,
    XtX_inv = XtX_inv,
    sigma2_hat = sigma2_hat,
    n = n,
    p = p
  )
}
pred_blm <- function(obj, Xt){
  fit <- obj
  Xtest <- Xt
  Xtest <- cbind(1, Xtest)
  beta_hat <- fit$beta_hat
  XtX_inv <- fit$XtX_inv
  sigma2_hat <- fit$sigma2_hat
  num_draws <- 1000
  predictions <- matrix(NA, nrow = num_draws, ncol = nrow(Xtest))
  for (i in 1:nrow(Xtest)) {
    mu_pred <- as.numeric(Xtest[i, ] %*% beta_hat)
    pred_var <- sigma2_hat * (1 + as.numeric(Xtest[i, ] %*% XtX_inv %*% Xtest[i, ]))
    predictions[, i] <- rnorm(num_draws, mean = mu_pred, sd = sqrt(pred_var))
  }
  predictions
}
model_names <- c(model_names, "blm")

# 5. Bayesian CART
library(tgp)
fit_bcart <- function(X, y){
  if(var(y) == 0) y <- y + rnorm(y, 0, 1e-6)
  out <- list(X=X, y=y)
  return(out)
}
pred_bcart <- function(obj, Xt){
  mod <- tgp::bcart(obj$X, obj$y, Xt, zcov=TRUE)
  yhat <- as.numeric(mod$ZZ.mean)           # mean vector (length n)
  Sigma <- mod$ZZ.s2                        # covariance matrix (n x n)
  n <- length(yhat)

  # Cholesky decomposition: Sigma = L %*% t(L)
  L <- chol(Sigma)

  # Generate nsamp samples of standard normals
  nsamp <- 1000
  z <- matrix(rnorm(nsamp * n), nrow = n)

  # Generate samples: Y = mu + L %*% z
  samples <- matrix(NA, nrow = nsamp, ncol = n)
  for(i in 1:nsamp){
    samples[i, ] <- yhat + L %*% z[,i]
  }

  return(samples)
}
model_names <- c(model_names, "bcart")

# 6. Bayesian Treed Linear Model
fit_btreelm <- function(X, y){
  if(var(y) == 0) y <- y + rnorm(y, 0, 1e-6)
  out <- list(X=X, y=y)
  return(out)
}
pred_btreelm <- function(obj, Xt){
  mod <- tgp::btlm(obj$X, obj$y, Xt, zcov=TRUE)
  yhat <- as.numeric(mod$ZZ.mean)           # mean vector (length n)
  Sigma <- mod$ZZ.s2                        # covariance matrix (n x n)
  n <- length(yhat)

  # Cholesky decomposition: Sigma = L %*% t(L)
  L <- chol(Sigma)

  # Generate nsamp samples of standard normals
  nsamp <- 1000
  z <- matrix(rnorm(nsamp * n), nrow = n)

  # Generate samples: Y = mu + L %*% z
  samples <- matrix(NA, nrow = nsamp, ncol = n)
  for(i in 1:nsamp){
    samples[i, ] <- yhat + L %*% z[,i]
  }

  return(samples)
}
model_names <- c(model_names, "btreelm")

# 7. Sparse Khaos
library(khaos)
library(glmnet)
fit_spce <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y,0,1e-6)
  sparse_khaos(X,y, verbose=FALSE)
}
pred_spce <- function(obj, Xt) predict(obj, Xt)
model_names <- c(model_names, "spce")

# 8. Adaptive Khaos
fit_apce <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y,0,1e-6)
  adaptive_khaos(X, y, verbose=FALSE)
}
pred_apce <- function(obj, Xt) predict(obj, Xt)
model_names <- c(model_names, "apce")

# 9. Conformal Random Forest
library(conforest)
fit_confrf <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y,0,1e-9)
  rfok(X, y)
}
pred_confrf <- function(obj, Xt) predict(obj, Xt)
model_names <- c(model_names, "confrf")

# 10. Sparse GP (Fully independent training conditionals; FITC)
library(gplite)
fit_fitcgp <- function(X, y){
  if(sd(y) == 0){
    y <- y + rnorm(y, 0, 1e-7)
  }else{
    y <- y + rnorm(y, 0, 1e-7*sd(y)) # Do this for robustness
  }
  n <- length(y)
  ni <- 4*floor(sqrt(n))
  ni <- max(5, ni)
  ni <- min(100, ni, round(n*0.75))
  gp <- gp_init(method=method_fitc(num_inducing = ni))
  gp <- gp_optim(gp, X, y, restarts=2)
  gp$y <- y
  return(gp)
}
pred_fitcgp <- function(obj, Xt){
  jitter <- 1e-6
  while(TRUE){
    t_pred <- try({
      gp_draw(obj, Xt, draws=1000, jitter=jitter)
    }, silent = TRUE)
    if(!inherits(t_pred, "try-error")){
      return(t(t_pred))
    }else{
      print(t_pred)
      jitter <- jitter * 10
      if(jitter >= 100*var(obj$y)) stop("gp draw failed even with jitter > var(y)")
    }
  }
}
model_names <- c(model_names, "fitcgp")

# 11. Basis GP (Random Fourier Features; Rahimi and Recht 2007)
fit_rffgp <- function(X, y){
  if(sd(y) == 0){
    y <- y + rnorm(y, 0, 1e-7)
  }else{
    y <- y + rnorm(y, 0, 1e-7*sd(y)) # Do this for robustness
  }
  n <- length(y)
  nb <- 2*floor(sqrt(n))
  nb <- min(512, nb)
  while(TRUE){
    gp <- try({
      gp1 <- gp_init(method=method_rf(num_basis=nb))
      gp1 <- gp_optim(gp1, X, y, restarts=2)
      gp1$y <- y
      gp1
    }, silent=TRUE)

    if(!inherits(gp, "try-error")){
      return(gp)
    }else{
      print(gp)
      if(nb <= 2) stop("numb basis 1 failed")
      nb <- floor(nb/2) + (floor(nb/2) %% 2)
    }
  }
}
pred_rffgp <- function(obj, Xt){
  jitter <- 1e-6
  while(TRUE){
    t_pred <- try({
      gp_draw(obj, Xt, draws=1000, jitter=jitter)
    }, silent = TRUE)
    if(!inherits(t_pred, "try-error")){
      return(t(t_pred))
    }else{
      print(t_pred)
      jitter <- jitter * 10
      if(jitter >= 100*var(obj$y)) stop("gp draw failed even with jitter > var(y)")
    }
  }
}
model_names <- c(model_names, "rffgp")

# 12. Local Approximate GP
library(laGP)
fit_lagp <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-6)
  if(max(table(y)) > 10) y <- y + rnorm(y, 0, 1e-6)
  list(X=X, y=y)
}
pred_lagp <- function(obj, Xt){
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
model_names <- c(model_names, "lagp")

# 13. Local Approximate GP (SoD)
fit_alcgp <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-6)
  list(X=X, y=y)
}
pred_alcgp <- function(obj, Xt){
  n <- nrow(Xt)
  nn <- max(100, min(300, n, 2 * round(sqrt(n))))
  mod <- laGP(Xref=Xt, obj$X, obj$y, start=10, end=nn)
  m <- mod$mean
  v <- mod$df
  s2 <- mod$s2 * (v-2) / v
  preds <- matrix(NA, nrow=1000, ncol=nrow(Xt))
  for(i in 1:nrow(Xt)){
    preds[,i] <- m[i] + rt(1000, v) * sqrt(s2[i])
  }
  return(preds)
}
model_names <- c(model_names, "alcgp")

# 14. Matching Pursuit SoD GP (Keerthi & Chu)
library(spareGParts)
library(GpGp)
library(GPvecchia) # These shouldn't be necesarry... but here we are.
library(GPfit)
library(cluster)
fit_mpgp <- function(X, y){
  if(sd(y) == 0){
    y <- y + rnorm(y, 0, 1e-7)
  }else{
    y <- y + rnorm(y, 0, sd(y)*1e-7)
  }

  n <- length(y)
  m <- min(1000, n-1, max(100, 2*floor(sqrt(n))))

  res <-  mpgp(X, y, m=m)
  return(res)
}
pred_mpgp <- function(obj, Xt){
  predict(obj, Xt)
}
model_names <- c(model_names, "mpgp")

# 15. Scaled Vecchia GP
fit_svecgp <- function(X, y){
  y_og <- y
  if (sd(y) == 0) y <- y + rnorm(y,0,1e-7)

  if(ncol(X) == 1){
    X <- cbind(X, rnorm(length(y), 1, 0.001)) # add a dummy column
  }

  nug <- NULL
  tries <- 0
  fit <- NULL

  while(TRUE){
    tries <- tries + 1
    fit <- try({
      svecgp(X, y, nug=nug)
    }, silent=TRUE)
    if(!inherits(fit, "try-error")){
      return(fit)
    }else{
      cat("svecgp failed with nug = ", nug, "on try", tries, "\n")
      cat("error says:", fit, "\n")
      if(is.null(nug)){
        nug <- var(y) / 100
      }else{
        nug <- nug * 10
        y <- y_og + rnorm(y, 0, sqrt(nug))
      }
    }
    if(tries >= 30) stop("too many nugget failures")
  }
}
pred_svecgp <- function(obj, Xt){
  if(ncol(Xt) == 1) Xt <- cbind(Xt, rnorm(nrow(Xt),1,0.001))
  predict(obj, Xt)
}
model_names <- c(model_names, "svecgp")

# 16. Bayesian NN
library(bnns)
fit_bnn <- function(X,  y){
  data <- as.data.frame(X)
  data$y <- y
  obj <- bnns(y~., data=data, L=2, nodes=c(8,8),
              iter=1200, warmup=200, thin=2, chains=2, cores=1)
  return(obj)
}
pred_bnn <- function(obj, Xt){
  preds <- t(predict(obj, as.data.frame(Xt)))
  return(preds)
}
model_names <- c(model_names, "bnn")

pred_bnn2 <- function(obj, Xt){
  preds <- t(predict(obj, as.data.frame(Xt)))
  preds <- rbind(preds, preds, preds)[1:1000,]
  return(preds)
}

# 17. T-BASS
library(GBASS)
fit_tbass <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-7)
  tbass(X, y, df=3)
}
pred_tbass <- function(obj, Xt) predict(obj, Xt)
model_names <- c(model_names, "tbass")

# 18. Median BASS
fit_qbass <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-7)
  qbass(X, y, q=0.5)
}
pred_qbass <- function(obj, Xt) predict(obj, Xt)
model_names <- c(model_names, "qbass")

# 19. Relevance Vector Machine
#library(spareGParts)
#library(glmnet)
fit_rvm <- function(X, y) rvm(X, y)
pred_rvm <- function(obj, Xt) predict(obj, Xt)
model_names <- c(model_names, "rvm")

# 20. NGBoost
library(ngboost)
use_condaenv("ngbenv", required = TRUE)
py_config()
fit_ngboost <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-7)
  model <- NGBRegression$new(Dist = Dist("Normal"),
                             Base = sklearner(),
                             Score = Scores("MLE"),
                             natural_gradient =TRUE,
                             n_estimators = 600,
                             learning_rate = 0.002,
                             minibatch_frac = 0.8,
                             col_sample = 0.9,
                             verbose = TRUE,
                             verbose_eval = 100,
                             tol = 1e-5)

  model$fit(X=X, Y=y)
  return(model)
}
pred_ngboost <- function(obj, Xt){
  samples <- 1000
  preds <- matrix(NA, nrow=samples, ncol=nrow(Xt))
  for(i in 1:nrow(Xt)){
    distt <- obj$pred_dist(Xt[i,,drop=FALSE])
    mu <- distt$mean()
    sigma <- distt$std()
    preds[,i] <- rnorm(samples, mu, sigma)
  }
  return(preds)
}
model_names <- c(model_names, "ngboost")

# 21. Bayesian LASSO
library(BayesianLasso)
fit_blasso <- function(X, y){
  # package always fails if ncol(X) <= 2
  # so add some dummy columns
  if(ncol(X) == 2){
    X <- cbind(rnorm(nrow(X), 0, 1), X)
  }
  if(ncol(X) == 1){
    X <- cbind(rnorm(nrow(X), 0, 1), X)
    X <- cbind(rnorm(nrow(X), 0, 1), X)
  }
  max_attempts = 5

  mu_y <- mean(y)
  sd_y <- sd(y)
  if(sd_y == 0){
    y <- y + rnorm(y, 0, 1e-6)
    mu_y <- mean(y)
    sd_y <- sd(y)
  }
  y_std <- (y - mu_y) / sd_y

  # Starting priors
  a1 <- 1
  b1 <- 1
  u1 <- 0.1
  v1 <- 0.1
  nsamples <- 3000
  beta_init <- rep(0.1, ncol(X))
  lambda_init <- 1
  sigma2_init <- 1

  attempt <- 1
  while (attempt <= max_attempts) {
    fit <- try(
      BayesianLasso::Modified_Hans_Gibbs(
        X, y_std, a1, b1, u1, v1,
        nsamples, beta_init, lambda_init, sigma2_init,
        verbose = FALSE
      ),
      silent = TRUE
    )
    if (!is.list(fit) || is.null(names(fit)) || !all(c("mBeta", "vsigma2") %in% names(fit))) {
      cat("blasso failure: class=", class(fit), "value=", capture.output(str(fit)), "\n", file=sprintf("fit_fail_%d.log", FUNC_IDX))
      saveRDS(list(fit=fit, X=X, y=y, TMPDIR=Sys.getenv("TMPDIR")), file=sprintf("fit_fail_%d.rds", FUNC_IDX))
      # fallback logic here
    }
    # Check if fit failed, or if all coefficients are zero/empty (super-sparse)
    fail <- inherits(fit, "try-error") ||
      is.null(fit$mBeta) || ncol(as.matrix(fit$mBeta)) == 0 ||
      all(colSums(abs(as.matrix(fit$mBeta))) == 0)
    if (!fail) {
      fit$mu_y <- mu_y
      fit$sd_y <- sd_y
      fit$a1 <- a1
      fit$b1 <- b1
      fit$attempts <- attempt
      fit$fallback <- FALSE
      class(fit) <- "blasso"
      return(fit)
    }
    # Update priors for next attempt
    a1 <- a1 / 2
    b1 <- b1 * 2
    lambda_init <- lambda_init / 4
    attempt <- attempt + 1
  }
  warning("Bayesian Lasso failed after multiple attempts. Falling back to mean-only model.")
  # Fallback: mean-only model
  fallback <- list(
    mBeta = matrix(0, nrow = nsamples, ncol = ncol(X)),
    vsigma2 = rep(var(y), nsamples),
    mu_y = mu_y,
    sd_y = sd_y,
    a1 = a1,
    b1 = b1,
    attempts = attempt - 1,
    fallback = TRUE
  )
  class(fallback) <- "blasso"
  return(fallback)
}
pred_blasso <- function(obj, Xt){
  # package always fails if ncol(X) <= 2
  # so add some dummy columns
  if(ncol(Xt) == 2){
    Xt <- cbind(rnorm(nrow(Xt), 0, 1), Xt)
  }
  if(ncol(Xt) == 1){
    Xt <- cbind(rnorm(nrow(Xt), 0, 1), Xt)
    Xt <- cbind(rnorm(nrow(Xt), 0, 1), Xt)
  }
  print("predicting now")
  nt <- nrow(Xt)
  ind <- seq(1001, 3000, by=2)
  preds <- matrix(NA, nrow=length(ind), ncol=nt)
  for(i in seq_along(ind)){
    beta <- obj$mBeta[i,]
    sigma <- sqrt(obj$vsigma2[i])
    preds[i,] <- rnorm(nt, Xt %*% beta, sigma)
  }
  preds <- obj$mu_y + obj$sd_y * preds
  return(preds)
}
model_names <- c(model_names, "blasso")


# 22. Full GP
library(hetGP)
fit_gp <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-7)
  mu_y <- mean(y)
  sigma_y <- sd(y)
  y <- (y-mu_y)/sigma_y

  fit <- hetGP::mleHomGP(X, y, noiseControl=list(g_max=1, g_min=1e-9))
  fit$mu_y <- mu_y
  fit$sigma_y <- sigma_y
  return(fit)
}
pred_gp <- function(obj, Xt){
  tmp <- predict(obj, Xt)
  mu <- tmp$mean
  sigma <- sqrt(abs(tmp$sd2 + tmp$nugs))
  preds <- matrix(NA, nrow=1000, ncol=nrow(Xt))
  for(i in 1:nrow(Xt)){
    preds[,i] <- rnorm(1000, mu[i], sigma[i])
  }
  preds <- obj$mu_y + obj$sigma_y * preds
}
model_names <- c(model_names, "gp")

# 23. Deep GP
library(deepgp)
fit_deepgp <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-7)
  mu_y <- mean(y)
  sigma_y <- sd(y)
  y <- (y-mu_y)/sigma_y
  cat("fitting")
  fit <- deepgp::fit_two_layer(X, y, nmcmc=1000, vecchia=TRUE)
  fit <- trim(fit, 500, 5)
  fit$mu_y <- mu_y
  fit$sigma_y <- sigma_y
  return(fit)
}

pred_deepgp <- function(obj, Xt){
  cat("predicting")
  tmp <- predict(obj, Xt)
  mu <- tmp$mean
  sigma <- sqrt(abs(tmp$s2))

  print("mu")
  print(mu)
  print("sigma")
  print(sigma)

  mu[is.na(mu)] <- obj$mu_y
  sigma[is.na(sigma)] <- obj$sigma_y

  preds <- matrix(NA, nrow=1000, ncol=nrow(Xt))
  for(i in 1:nrow(Xt)){
    preds[,i] <- rnorm(1000, mu[i], sigma[i])
  }
  preds <- obj$mu_y + obj$sigma_y * preds
}
model_names <- c(model_names, "deepgp")

# 24. Tree GP
library(tgp)
fit_treegp <- function(X, y){
  if(var(y) == 0) y <- y + rnorm(y, 0, 1e-6)
  out <- list(X = X, y = y)
  return(out)
}

pred_treegp <- function(obj, Xt){
  mod <- tgp::btgp(obj$X, obj$y, Xt, zcov = TRUE)
  yhat <- as.numeric(mod$ZZ.mean)   # mean vector (length n)
  Sigma <- mod$ZZ.s2                # covariance matrix (n x n)
  n <- length(yhat)
  nsamp <- 1000

  # Cholesky decomposition
  nug_vec <- c(0, 1e-7, 1e-4, 1e-1)
  for(i in seq_along(nug_vec)){
    nugI <- diag(rep(nug_vec[i], n))
    L <- try(chol(Sigma + nugI), silent=TRUE)
    if(!inherits(L, "try-error")){
      break
    }
  }

  # Sample with cholesky
  z <- matrix(rnorm(nsamp * n), nrow = n)
  samples <- t(yhat + crossprod(L, z))
  return(samples)
}
model_names <- c(model_names, "treegp")

# 25. Het GP
library(hetGP)
fit_hetgp <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-7)
  mu_y <- mean(y)
  sigma_y <- sd(y)
  y <- (y-mu_y)/sigma_y

  fit <- hetGP::mleHetGP(X, y, noiseControl=list(g_max=1, g_min=1e-9))
  fit$mu_y <- mu_y
  fit$sigma_y <- sigma_y
  return(fit)
}
pred_hetgp <- function(obj, Xt){
  tmp <- predict(obj, Xt)
  mu <- tmp$mean
  sigma <- sqrt(abs(tmp$sd2 + tmp$nugs))
  preds <- matrix(NA, nrow=1000, ncol=nrow(Xt))
  for(i in 1:nrow(Xt)){
    preds[,i] <- rnorm(1000, mu[i], sigma[i])
  }
  preds <- obj$mu_y + obj$sigma_y * preds
}
model_names <- c(model_names, "hetgp")

# 26. Random Forest + Bagging (bootstrap)
library(randomForest)
fit_bootrf <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-7)
  n_bags = 100
  bag_models <- vector("list", n_bags)
  n <- nrow(X)
  for (b in 1:n_bags) {
    idx <- sample(n, n, replace = TRUE)
    bag_models[[b]] <- randomForest::randomForest(
      x = X[idx, , drop = FALSE],
      y = y[idx])
  }
  bag_models
}

pred_bootrf <- function(obj, Xt){
  n_bags <- length(obj)
  n_test <- nrow(Xt)
  preds_mat <- matrix(NA, nrow = n_bags, ncol = n_test)
  for (b in 1:n_bags) {
    preds_mat[b, ] <- predict(obj[[b]], newdata = Xt)
  }
  preds_mat
}
model_names <- c(model_names, "bootrf")

# 27. Bayesian Committee Machine
library(spareGParts)
fit_bcmgp <- function(X, y){
  if(sd(y) == 0){
    y <- y + rnorm(y, 0, 1e-7)
  }else{
    y <- y + rnorm(y, 0, sd(y)*1e-7)
  }
  bcmgp(X, y)
}
pred_bcmgp <- function(obj, Xt){
  predict(obj, Xt)
}
model_names <- c(model_names, "bcmgp")

# 28. Robust GaSP
library(RobustGaSP)
fit_rgasp <- function(X, y){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-7)
  print("trying gasp")
  RobustGaSP::rgasp(X, y, nugget.est=TRUE)
}
pred_rgasp <- function(obj, Xt){
  print("predicitng gasp")
  out <- predict.rgasp(obj, Xt)
  mu <- out$mean
  s  <- out$sd
  preds <- matrix(NA, nrow=1000, ncol=length(mu))
  for(i in seq_along(mu)){
    preds[,i] <- rnorm(1000, mu[i], s[i])
  }
  print("returning")
  return(preds)
}
model_names <- c(model_names, "rgasp")

# 29. Baseline (fallback model)
fit_baseline <- function(X, y){
  stop("Use fallback model")
}
pred_baseline <- function(obj, Xt){
  stop("Shouldn't even get here")
}
model_names <- c(model_names, "baseline")
