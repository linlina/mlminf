source("mlminf.R")

# Load in data.
ribo <-
  data.frame(read.csv(
    "Data/riboflavingrouped.csv",
    header = TRUE,
    row.names = 1
  ))
ribo_group <-
  read.csv("Data/riboflavingrouped_structure.csv", header = FALSE)

# Clean up and process data.
y <- as.numeric(ribo[1,])
riboX <- ribo[-1, ]

cell_labels <- names(riboX)
genes <- rownames(riboX)

X <- as.data.frame(t(riboX))
rownames(X) <- cell_labels
names(X) <- genes
N <- nrow(X)

group_id <- as.numeric(ribo_group[, 1])
X_sd <- scale(X)
y_sd <- scale(y)
cors <- crossprod(X_sd, y_sd)

X_select <- X_sd
genes_select <- genes

# Set up inputs to method functions.
m <- length(unique(group_id))
n_list <- rep(0, m)
group_tab <- as.numeric(table(group_id))
for (i in 1:m) {
  n_list[i] <- group_tab[unique(group_id)[i]]
}

N <- sum(n_list)
p <- ncol(X_select)
q <- 1

Z_list <- list()
for (i in 1:m) {
  Z_list[[i]] <- rep(1, n_list[i])
}
Z_block <- bdiag(Z_list)

# Run our method.
ribo_mlm <-
  MLMInf(
    N,
    m,
    n_list,
    p,
    q,
    as.matrix(X_select),
    Z_block,
    y_sd,
    lambda_ridge = 1 / N,
    known = FALSE,
    sims = 10000
  )

# Assuming fully independent structure.
ribo_lm <-
  LMInf(N,
        p,
        as.matrix(X_select),
        y_sd,
        lambda_ridge = 1 / N,
        sims = 10000)
