source("mlminf.R")
library(corpcor)

GenerateData <- function(N, n, m, p, q, X, beta, Z, v, tau, seed) {
  set.seed(seed)
  y <- rep(0, N)
  resid <- rep(0, N)
  epsilon <- rnorm(N, sd = v)
  alpha <- matrix(0, m, q)
  for (i in 1:m) {
    alpha[i,] <- rnorm(q, mean = 0, sd = tau)
    resid[((i - 1) * n + 1):(i * n)] <-
      tcrossprod(Z[((i - 1) * n + 1):(i * n),], t(alpha[i, ])) + epsilon[((i -
                                                                             1) * n + 1):(i * n)]
    y[((i - 1) * n + 1):(i * n)] <-
      tcrossprod(X[((i - 1) * n + 1):(i * n),], t(beta)) + resid[((i - 1) * n +
                                                                    1):(i * n)]
  }
  r <- list(y = y, alpha = alpha, resid = resid)
  return(r)
}


# To generate data for Figures 1 and 2 (for plot code, see process.R).
m <- 25
n <- 6
N <- m * n

q_list <- c(1, 2)
s_list <- c(5, 10)
beta_min_list <- c(0.5, 1)
model_list <- 1:2

p_list <- c(300, 600)

for (p_index in 1:length(p_list)) {
  for (q_index in 1:length(q_list)) {
    for (s_index in 1:length(s_list)) {
      for (beta_min_index in 1:length(beta_min_list)) {
        for (model_index in 1:length(model_list)) {
          p <- p_list[p_index]
          q <- q_list[q_index]
          s <- s_list[s_index]
          beta_min <- beta_min_list[beta_min_index]
          model <- model_list[model_index]
          
          
          if (model == 1) {
            rho <- 0.2
            times <- 1:(p + q)
            H <- abs(outer(times, times, "-"))
            Sigma <- rho ^ H
          } else {
            Sigma <- diag(1, p + q)
          }
          
          set.seed(1)
          Data <- rmvnorm(N, mean = rep(0, p + q), sigma = Sigma)
          X <- Data[, 1:p]
          X <- sweep(X, 2, sqrt(colSums(X ^ 2)), '/') * sqrt(N)

          Z <- Data[, -(1:p), drop = FALSE]
          
          Z_block  <- matrix(0, N, m * q)
          for (i in 1:m) {
            Z_block[((i - 1) * n + 1):(i * n), ((i - 1) * q + 1):(i * q)] <-
              Z[((i - 1) * n + 1):(i * n), ]
          }
          
          Z_block <-
            sweep(Z_block, 2, sqrt(colSums(Z_block ^ 2)), '/') * sqrt(n)
          
          Z <- matrix(0, N, q)
          for (i in 1:m) {
            Z[((i - 1) * n + 1):(i * n), ] <-
              Z_block[((i - 1) * n + 1):(i * n), ((i - 1) * q + 1):(i * q)]
          }
          
          v <- 0.5
          tau <- 1
          
          V <- v ^ 2 * diag(1, N) + tau ^ 2 * tcrossprod(Z_block)
          maxV <- max(eigen(V)$values)
          lambda0 = sqrt(2 * log(p) / N * maxV)
          
          trials <- 1000
          
          beta <- c(rep(beta_min, s), rep(0, p - s))
          test_tracker <- pval_tracker <- matrix(0, trials, p)
          
          seed <- 0
          for (index in 1:trials) {
            seed <- seed + 1
            data <-
              GenerateData(N, n, m, p, q, X, beta, Z, v, tau, seed)
            y <- data$y
            run <-
              tryCatch({
                MLMInf(
                  N,
                  m,
                  n,
                  p,
                  q,
                  X,
                  Z,
                  y,
                  lambda_ridge = 1 / N,
                  known = FALSE,
                  sims = 10000,
                  lambda0 = lambda0
                )
              }, error = function(err) {
                return(1)
              })
            while (length(run) == 1) {
              seed <- seed + 1
              data <-
                GenerateData(N, n, m, p, q, X, beta, Z, v, tau, seed)
              y <- data$y
              run <-
                tryCatch({
                  MLMInf(
                    N,
                    m,
                    n,
                    p,
                    q,
                    X,
                    Z,
                    y,
                    lambda_ridge = 1 / N,
                    known = FALSE,
                    sims = 10000,
                    lambda0 = lambda0
                  )
                }, error = function(err) {
                  return(1)
                })
            }
            pval_tracker[index,] <- run
            test_tracker[index,] <- (run <= 0.05) * 1
          }
          
          filename <-
            paste(
              paste(
                "Output/Ridge_lmm",
                "p",
                p,
                "q",
                q,
                "model",
                model,
                "s",
                s,
                "beta.min",
                beta.min,
                "fixed_adjust_delta_multi",
                sep = "_"
              ),
              ".RData",
              sep = ""
            )
          save(p,
               q,
               model,
               s,
               beta_min,
               test_tracker,
               tau_tracker,
               v_tracker,
               file = filename)
          
        }
      }
    }
  }
}

# To generate data for Figure 3 (for plot code, see process.R).
model <- 1
p <- 300
q <- 1

m <- 25
n <- 6
N <- m * n

rho <- 0.2
times <- 1:(p + q)
H <- abs(outer(times, times, "-"))
Sigma <- rho ^ H


# Generate data.
set.seed(1)
Data <- rmvnorm(N, mean = rep(0, p + q), sigma = Sigma)
X <- Data[, 1:p]
X <- sweep(X, 2, sqrt(colSums(X ^ 2)), '/') * sqrt(N)
pX <- pseudoinverse(crossprod(X))

Z <- Data[, -(1:p), drop = FALSE]

Z_block  <- matrix(0, N, m * q)
for (i in 1:m) {
  Z_block[((i - 1) * n + 1):(i * n), ((i - 1) * q + 1):(i * q)] <-
    Z[((i - 1) * n + 1):(i * n), ]
}

Z_block <-  sweep(Z_block, 2, sqrt(colSums(Z_block ^ 2)), '/') * sqrt(n)

Z <- matrix(0, N, q)
for (i in 1:m) {
  Z[((i - 1) * n + 1):(i * n), ] <-
    Z_block[((i - 1) * n + 1):(i * n), ((i - 1) * q + 1):(i * q)]
}

v <- 0.5
tau <- 1

V <- v ^ 2 * diag(1, N) + tau ^ 2 * tcrossprod(Z_block)
maxV <- max(eigen(V)$values)
lambda0 = sqrt(2 * log(p) / N * maxV)

trials <- 1000
s <- 5
beta <- c(0.05, 2, 4, 3, 0.1, rep(0, p - s))

lmm_tracker <- rep(0, p)
for (index in 1:trials) {
  data <- GenerateData(N, n, m, p, q, X, beta, Z, v, tau, index)
  y <- data$y
  run <-
    MLMInf(
      N,
      m,
      n,
      p,
      q,
      X,
      Z,
      y,
      lambda_ridge = 1 / N,
      known = FALSE,
      sims = 10000,
      lambda0 = lambda0
    )
  temp <- (beta >= run$l & beta <= run$u) * 1
  lmm_tracker <- lmm_tracker + temp / trials
}

ridge_tracker <- rep(0, p)
trials <- 1000 
for (index in 1:trials) {
  data <- GenerateData(N, n, m, p, q, X, beta, Z, v, tau, index)
  y <- data$y
  run <- LMInf(N, p, X, y, lambda_ridge = 1 / N, sims = 10000)
  temp <- (beta >= run$l & beta <= run$u) * 1
  ridge_tracker <- ridge_tracker + temp / trials
}

lasso_tracker <- rep(0, p)
lambda_list <- exp(seq(-4, 0, length = 100))
Theta_hat <-
  MakeThetaHat(N, p, X, lambda_list, K = 10, type = "1se")
for (index in 1:trials) {
  data <- GenerateData(N, n, m, p, q, X, beta, Z, v, tau, index)
  y <- data$y
  run <-
    LassoLM(N,
            p,
            X,
            y,
            Theta_hat = Theta_hat,
            sims = 10000,
            type = "1se")
  temp <- (beta >= run$l & beta <= run$u) * 1
  lasso_tracker <- lasso_tracker + temp / trials
}

filename <-
  paste(
    paste("Output/Ridge_lmm", "conf_int", "fixed", "ridgedelta", sep = "_"),
    ".RData",
    sep = ""
  )
save.image(file = filename)
