library(devtools)
load_all("..")

beta_pmf <- function(v, theta) {
  probs <- dbeta(v, shape1 = theta$mu * theta$lambda, shape2 = (1 - theta$mu) * theta$lambda)
  probs / sum(probs)
}

initialize_theta_beta <- function(C, v, K, hparams) {
  vmin <- min(v)
  vmax <- max(v)
  mus <- seq(vmin, vmax, length.out = K + 2)[2:(K + 1)]
  lapply(seq_len(K), function(k) {
    list(mu = mus[k], lambda = 25)
  })
}

update_theta_beta <- function(C, v, Gamma, theta, hparams) {
  K <- dim(Gamma)[3]
  lapply(seq_len(K), function(k) {
    weights <- colSums(C * Gamma[, , k, drop = FALSE][, , 1])
    total <- sum(weights)
    if (total <= 0) {
      return(theta[[k]])
    }

    mu_hat <- sum(weights * v) / total
    var_hat <- sum(weights * (v - mu_hat)^2) / total
    lambda_hat <- if (var_hat <= 0) {
      hparams$lambda_bounds[2]
    } else {
      mu_hat * (1 - mu_hat) / var_hat - 1
    }
    lambda_hat <- max(hparams$lambda_bounds[1], min(hparams$lambda_bounds[2], lambda_hat))
    mu_hat <- max(hparams$mu_eps, min(1 - hparams$mu_eps, mu_hat))

    list(mu = mu_hat, lambda = lambda_hat)
  })
}

set.seed(1)
v <- seq(0.05, 0.95, by = 0.1)
true_theta <- list(
  list(mu = 0.25, lambda = 18),
  list(mu = 0.75, lambda = 30)
)
true_w <- rbind(
  c(0.85, 0.15),
  c(0.55, 0.45),
  c(0.20, 0.80)
)

pmf_mat <- do.call(rbind, lapply(true_theta, function(th) beta_pmf(v, th)))
counts_total <- 1500
C <- t(apply(true_w %*% pmf_mat, 1, function(p) as.vector(rmultinom(1, counts_total, p))))

fit <- pmfmix(
  C = C,
  v = v,
  K = 2,
  f = beta_pmf,
  initialize_theta = initialize_theta_beta,
  update_theta = update_theta_beta,
  hparams = list(alpha = c(1, 1), mu_eps = 1e-6, lambda_bounds = c(2, 200)),
  control = list(nstart = 3, niter = 25, abstol = 1e-6),
  fixed = list(w = NULL, theta = NULL, Gamma = NULL),
  verbose = TRUE
)

cat("Mixture weights:\n")
print(round(fit$params$w, 3))
cat("Component parameters:\n")
print(lapply(fit$params$theta, function(th) list(mu = round(th$mu, 3), lambda = round(th$lambda, 3))))
