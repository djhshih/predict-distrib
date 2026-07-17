library(devtools) 
load_all("..")

beta_lpmf <- function(v, theta) {
  dbeta(v, shape1 = theta$mu * theta$lambda, shape2 = (1 - theta$mu) * theta$lambda, log=TRUE)
}

initialize_theta_beta <- function(C, v, K, hparams) {
  vmin <- min(v)
  vmax <- max(v)
  mus <- seq(vmin, vmax, length.out = K + 2)[2:(K + 1)]
  lapply(seq_len(K), function(k) {
    list(mu = mus[k], lambda = 1)
  })
}

update_theta_beta <- function(C, v, params, hparams) {
  N <- nrow(C);
  J <- ncol(C);
  K <- ncol(params$W);

  # 2*K parameters in theta to optimize (mu and lambda)
  mparam_transform <- function(a) {
    mu <- logistic(a[1:K]);
    lambda <- exp(a[(K+1):(K+K)]);
    list(mu = mu, lambda = lambda)
  }
  
  mparam_rev_transform <- function(theta) {
    mu <- unlist(lapply(theta, function(th) th$mu));
    lambda <- unlist(lapply(theta, function(th) th$lambda));
    c(logit(mu), log(lambda))
  }

  # Transform activities a to a probability mass function that is evaluated at xs
  # return N x M matrix, where each row is a probability mass function
  lwf_transform <- function(a) {
    theta <- mparam_transform(a);
    lp <- unlist(lapply(v,
      # mixture of beta distributions
      function(x) log(t(params$W)) + beta_lpmf(x, theta)
    ));
    # W^T is K by N,  dbeta(x, ...) is K   ->  each item is K by N
    # output is K by N by J; need N by J by K
    aperm(array(lp, c(K, N, J)), c(2, 3, 1))
  }

  # - E_Z[ log p(C, Z, theta) ]
  objective_q <- function(a) {
    - sum( params$Z * lwf_transform(a) )
  }

  lpdf_transform <- function(a) {
    theta <- mparam_transform(a);
    lp <- unlist(lapply(v,
      # mixture of beta distributions
      # numeric overflow can occur due to small dbeta, causing log(0) = -Inf
      function(x) apply(t(log(params$W)) + beta_lpmf(x, theta), 2, matrixStats::logSumExp)
    ));
    # output is N by K
    matrix(lp, nrow=nrow(C))
  }

  # - log p(C, W, theta)
  objective_marginal <- function(a) {
    - sum( C * lpdf_transform(a) )
  }

  objective <- objective_q;
  #objective <- objective_marginal;

  a0 <- mparam_rev_transform(params$theta);
  # NB bounds are adjusted to avoid numeric overflow (-Inf and Inf)
  opt <- optim(a0, objective, method="L-BFGS-B", lower=-10, upper=10);
  theta <- mparam_transform(opt$par);

  lapply(seq_len(K), function(k) {
    list(mu = theta$mu[k], lambda = theta$lambda[k])
  })
}

set.seed(1)
v <- seq(0.05, 0.95, by = 0.05)
true_theta <- list(
  list(mu = 0.25, lambda = 10),
  list(mu = 0.75, lambda = 30)
)
true_w <- rbind(
  c(0.85, 0.15),
  c(0.55, 0.45),
  c(0.20, 0.80),
  c(0.40, 0.60),
  c(0.10, 0.90),
  c(0.90, 0.10)
)

pmf_mat <- do.call(rbind, lapply(true_theta, function(th) exp(beta_lpmf(v, th))))
pmf_mat <- pmf_mat / rowSums(pmf_mat);

target <- true_w %*% pmf_mat;
rowSums(target)
counts_total <- 1e5;
C <- t(apply(target, 1, function(p) as.vector(rmultinom(1, counts_total, p))))

fit <- pmfmix(
  C = C,
  v = v,
  K = 2,
  lf = beta_lpmf,
  initialize_theta = initialize_theta_beta,
  update_theta = update_theta_beta,
  hparams = list(alpha = c(1, 1)),
  control = list(nstart = 5, niter = 25, abstol = 1e-6),
  fixed = list(W = NULL, theta = NULL, Gamma = NULL),
  verbose = TRUE
)

plot(target, C / rowSums(C))

pmf_mat_hat <- exp(fit$params$lF);
plot(pmf_mat, pmf_mat_hat)
cor(c(pmf_mat), c(pmf_mat_hat))

target.hat <- fit$params$W %*% pmf_mat_hat;

print(target)
print(target.hat)

plot(target, target.hat)
cor(t(target), t(target.hat))
cor(c(target), c(target.hat))

par(mfrow=c(3, 2))
for (i in 1:nrow(target)) {
  plot(v, target[i, ], type="l")
  lines(v, target.hat[i, ], col="blue")
}

cat("Mixture weights:\n")
print(round(fit$params$W, 3))
cat("Component parameters:\n")
print(lapply(fit$params$theta, function(th) list(mu = round(th$mu, 3), lambda = round(th$lambda, 3))))

