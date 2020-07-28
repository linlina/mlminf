# Load libraries.
library(Matrix)
library(glmnet)
library(mvtnorm)
library(emulator)
library(MASS)
library(scalreg)


# Our method.
mlminf <- function(N, m, n_list, p, q, X, Z_block, y, lambda_ridge = 1/N, known = TRUE, sims = 10000, ...){
  svdX <- svd(X)
  pX <- tcrossprod(svdX$v)
  pX0 <- pX
  diag(pX0) <- 0
  
  eZ <- eigen(tcrossprod(Z_block))
  LZ <- eZ$vectors
  
  # Heuristic method for adjusting penalty parameter for Lasso.
  eZ_val <- eZ$values[which(abs(eZ$values) > 1e-5)]
  adjust <- sqrt(max(eZ_val)/mean(eZ_val))
  
  Sigma_hat <- crossprod(X)/N
  eX <- eigen(Sigma_hat)
  
  beta_ridge <- solve(Sigma_hat + 1/N*diag(1, p)) %*% crossprod(X, y)/N
  # If optimal Lasso penalty parameter value is known (unlikely).
  if (known) {
    lasso <- glmnet(X, y, family = "gaussian", lambda = lambda0, standardize = FALSE, intercept = FALSE)
  } else {
    sc <- scalreg(X, y)
    lambda_lasso <- adjust*sqrt(2*log(p)/N)*sc$hsigma
    lasso <- glmnet(X, y, family = "gaussian", lambda = lambda_lasso, standardize = FALSE, intercept = FALSE)
  }
  
  # Guessing the "support" of beta.
  active_lasso <- which(as.vector(lasso$beta) != 0)

  if (length(active.lasso) > 0){
    X1 <- X[, active.lasso]
    U <- cbind(X1, as.matrix(Z_block))
    R1 <- quad.tform.inv(crossprod(U), t(y) %*% U)
    R2 <- quad.tform.inv(crossprod(U), t(y) %*% U) - quad.tform.inv(crossprod(X1), t(y)%*%X1)
    v_hat <- sqrt((sum(y^2) - R1)/(N - ncol(U)))
    denom <- sum(diag(quad.form(diag(1, N) - quad.tform.inv(crossprod(X1), X1), as.matrix(Z_block))))
    tau_hat <- sqrt((R2 - v.hat^2 * ncol(Z_block))/denom)
  } else {
    R1 <- quad.tform.inv(crossprod(as.matrix(Z_block)), as.matrix(Z_block))
    v_hat <- sqrt(1/(N - sum(diag(R1)))*(as.numeric(crossprod(y)) - as.numeric(quad.form(R1, y))))
    tau_hat <- sqrt((crossprod(y) - N * v_hat^2)/N)
  }
  

  
  D <- diag(1/sqrt(v.hat^2 + tau.hat^2 * eZ$values))
  D2 <- diag(sqrt(v.hat^2 + tau.hat^2 * eZ$values))
  
  Ly <- quad.tform(D, LZ) %*% y 
  LX <- quad.tform(D, LZ) %*% X
  
  if (length(active.lasso) > 0){
    test <- lm(y ~ X[, active_lasso] - 1)
    beta.init <- rep(0, p)
    beta.init[active.lasso] <- test$coefficients
  } else {
    beta.init <- rep(0, p)
  }
  
  beta.corr <- (beta.ridge - crossprod(pX0, beta.init))
  
  Usims <- matrix(0, sims, p)
  Usims_p <- rep(0, sims)

  A1 <- quad.tform(D2, LZ)
  A2 <-  1/N*tcrossprod(quad.tform(diag(1/(eX$values + lambda_ridge)), eX$vectors), X)
  Delta <- apply(abs(pX0), 2, max)/diag(pX)*(log(p)/N)^(1/2 - 0.005)
  
  diagu <- sqrt(diag(tcrossprod(A2 %*% A1)))
  
  for (i in 1:sims) {
    epsilon <- rnorm(N, 0, 1)
    Lepsilon <- A1 %*% epsilon
    U <- A2 %*% Lepsilon
    Usims[i, ] <- U/diag(pX)
    Usims_p[i] <- min(2 * (1 - pnorm(as.vector(abs(U)/diagu))))
  }
  
  pval <- 2 * (1 - pnorm(pmax((abs(beta_corr/diagu) - Delta/diagu), 0)))
  sims_ecdf <- ecdf(Usims_p)
  pval_corr <- sims_ecdf(pval)
  
  lower <- apply(Usims, 2, quantile, prob = 0.025)
  upper <- apply(Usims, 2, quantile, prob = 0.975)
  
  return(list(pval = pval, pval_corr = pval_corr,
              l = beta.corr + lower - Delta, u = beta.corr + upper + Delta))
}

# Assuming independence, ridge method for linear fixed effect model.
lminf <- function(N, p, X, y, lambda_ridge = 1/N, sims = 10000){
  svdX <- svd(X)
  pX <- tcrossprod(svdX$v)
  pX0 <- pX
  diag(pX0) <- 0
  
  Sigma_hat <- crossprod(X)/N
  eX <- eigen(Sigma_hat)
  
  beta_ridge <- solve(Sigma_hat + 1/N*diag(1, p)) %*% crossprod(X, y)/N
  sc <- scalreg(X, y)
  
  beta_init <- sc$coefficients
  sigma_hat <- sc$hsigma
  
  beta_corr <- (beta_ridge - crossprod(pX0, beta_init))
  
  Usims <- matrix(0, sims, p)
  Usims_p <- rep(0, sims)
  A <- 1/N*tcrossprod(solve(Sigma_hat + lambda_ridge*diag(1,p)), X)
  Delta <- apply(abs(pX0), 2, max)/diag(pX)*(log(p)/N)^(1/2 - 0.005)
  
  diagu <- sqrt(diag(tcrossprod(A)))
  
  for (i in 1:sims){
    epsilon <- rnorm(N, 0, sigma_hat)
    U <- A %*% epsilon
    Usims[i, ] <- U/diag(pX)
    Usims_p[i] <- min(2 * (1 - pnorm(as.vector(abs(U)/diagu))))
  }
  
  pval <- 2 * (1 - pnorm(pmax((abs(beta_corr/diagu) - Delta/diagu),0)))
  temp <- ecdf(Usims_p)
  pval_corr <- temp(pval)
  
  lower <- apply(Usims, 2, quantile, prob = 0.025)
  upper <- apply(Usims, 2, quantile, prob = 0.975)
  
  return(list(pval = pval, pval_corr = pval_corr,
              l = beta.corr + lower - Delta, u = beta.corr + upper + Delta))
}

# Assuming independence, ridge method for linear fixed effect model.
lassolm <- function(N, p, X, y, Theta_hat, sims = 10000, type = c("min", "1se")){
  
  cvlasso <- cv.glmnet(X, y, family = "gaussian", nlambda = 100, nfolds = 10)
  
  if (type == "min"){
    best_lambda <- cvlasso$lambda.min
  } 
  if (type == "1se"){
    best_lambda <- cvlasso$lambda.1se
  }
  
  lasso <- glmnet(X, y, family = "gaussian", lambda = best.lambda)
  
  beta_hat <- as.vector(lasso$beta)
  sigma_hat <- sqrt(sum((y - X%*%beta_hat)^2)/(N - length(which(beta_hat != 0))))
  k_hat <- crossprod(X, y - X %*% lasso$beta)/(N*best_lambda)
  
  beta_corr <- beta_hat + best_lambda*Theta_hat%*%k_hat
  
  A <- 1/N*tcrossprod(Theta_hat, X)
  
  Usims <- matrix(0, sims, p)
  for (i in 1:sims){
    epsilon <- rnorm(N, 0, sigma_hat)
    U <- A %*% epsilon
    Usims[i,] <- U
  }
  
  lower <- apply(Usims, 2, quantile, prob = 0.025)
  upper <- apply(Usims, 2, quantile, prob = 0.975)
  
  return(list(l = beta_corr + lower, u = beta_corr + upper))
}

Theta_lambda_unit <- function(N, coord, X, data_select, lambda_list, K){
  error <- matrix(0, length(lambda_list), K)
  for (j in 1:K){
    whichj <- data_select == j
    temp <- glmnet(X[!whichj, -coord, drop=FALSE], X[!whichj, coord, drop=FALSE], lambda=lambda_list)
    predictions <- predict(temp, X[whichj, -coord, drop = FALSE], s = lambda_list)
    error[, j] <- apply((X[whichj, coord] - predictions)^2, 2, mean)
  }
  return(error)
}

Theta_lambda <- function(N, p, X, lambda_list, K = 10){
  data_select <- sample(rep(1:K, length = N))
  error <- mapply(Theta_lambda_unit,
                  coord = 1:p,
                  K = K,
                  data_select = list(data_select = data_select),
                  X = list(X = X),
                  lambda_list = list(lambda_list = lambda_list))
  err_array  <- array(unlist(error), dim = c(length(lambda_list), K, p))
  err_mean   <- apply(err_array, 1, mean)
  err_se     <- apply(apply(err_array, c(1, 2), mean), 1, sd)/sqrt(K)
  
  
  lambda_min <- lambda_list[which.min(err_mean)]
  lambda_min_sd <- err_se[which.min(err_mean)]
  
  lambda_1se <- max(lambda_list[err_mean < (min(err_mean) + lambda_min_sd)])
  return(list(lambda_min = lambda_min, lambda_1se = lambda_1se))
}

make_Theta_hat <- function(N, p, X, lambda_list, K = 10, type = c("min", "1se")){
  Beta_hat <- matrix(0, p, p)
  tau_hat_sq <- rep(0, p)
  
  if (type == "min"){
    best_lambda <- Theta_lambda(N, p, X, lambda_list, K = 10)$lambda_min
  }
  if (type == "1se"){
    best_lambda <- Theta_lambda(N, p, X, lambda_list, K = 10)$lambda_1se
  }
  
  for (j in 1:p){
    temp <- glmnet(X[, -j, drop = FALSE], X[, j, drop = FALSE], lambda = best_lambda)
    Beta_hat[-j,j] <- as.vector(temp$beta)
    tau_hat_sq[j] <- sum((X[,j] - X[,-j] %*% Beta_hat[j,-j])^2)/N + best_lambda*sum(abs(Beta_hat[j,-j]))
  }
  diag(Beta_hat) <- 1
  Theta_hat <- diag(1/tau_hat_sq) %*% Beta_hat 
  
  return(Theta_hat)
}

