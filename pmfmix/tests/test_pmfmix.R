library(devtools)
load_all("..")

beta_lpmf <- function(v, theta) {
  dbeta(v, shape1 = theta$mu * theta$lambda, shape2 = (1 - theta$mu) * theta$lambda, log=TRUE)
}

initialize_theta_beta <- function(C, v, K, hparams) {
  mu <- seq(0.2, 0.8, length.out = K)
  list(mu = mu, lambda = rep(K, 20))
}

update_theta_identity <- function(C, v, params, hparams) theta

assert <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

set.seed(42)
v <- seq(0.05, 0.95, by = 0.1)
K <- 2;
theta <- list(mu = c(0.3, 0.7), lambda = c(20, 25));
true_w <- rbind(c(0.8, 0.2), c(0.25, 0.75))
pmf_mat <- matrix(unlist(lapply(v, function(x) exp(beta_lpmf(x, theta)))), nrow=K);
C <- t(apply(true_w %*% pmf_mat, 1, function(p) as.vector(rmultinom(1, 4000, p))))

fit <- pmfmix(
  C = C,
  v = v,
  K = K,
  lf = beta_lpmf,
  initialize_theta = initialize_theta_beta,
  update_theta = update_theta_identity,
  control = list(nstart = 1, niter = 50, abstol = 1e-8),
  verbose = FALSE
)

assert(inherits(fit, "pmfmix"), "fit must inherit class 'pmfmix'")
assert(all(dim(fit$params$W) == c(2, 2)), "w must be a 2x2 matrix")
assert(max(abs(rowSums(fit$params$W) - 1)) < 1e-8, "rows of w must sum to 1")
assert(fit$params$W[1, 1] > fit$params$W[1, 2], "sample 1 should favor component 1")
assert(fit$params$W[2, 2] > fit$params$W[2, 1], "sample 2 should favor component 2")

cat("pmfmix tests passed\n")
