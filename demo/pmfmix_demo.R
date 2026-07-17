library(devtools) 
load_all("..")

beta_pmf <- function(v, theta) {
  dbeta(v, shape1 = theta$mu * theta$lambda, shape2 = (1 - theta$mu) * theta$lambda)
}

initialize_theta_beta <- function(C, v, K, hparams) {
  vmin <- min(v)
  vmax <- max(v)
  mus <- seq(vmin, vmax, length.out = K + 2)[2:(K + 1)]
  lapply(seq_len(K), function(k) {
    list(mu = mus[k], lambda = 25)
  })
}

logistic <- function(x) {
	1 / (1 + exp(-x))
}

logit <- function(x) {
  log(x) - log(1 - x)
}

update_theta_beta <- function(C, v, params, hparams) {
  N <- nrow(C);
  J <- ncol(C);
  K <- ncol(params$w);

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
  lpdf_transform <- function(a) {
    theta <- mparam_transform(a);
    lp <- with(theta, unlist(lapply(v,
      # mixture of beta distributions
      function(x) log(t(params$w) * dbeta(x, mu*lambda, (1 - mu)*lambda))
    )));
    # w^T is K by N,  dbeta(x, ...) is K   ->  each item is K by N
    # output is K by N by J; need N by J by K
    lp <- aperm(array(lp, c(K, N, J)), c(2, 3, 1))
    # message("lp:")
    # print(str(exp(lp) / sum(exp(lp))))
    lp
  }

  # negative log likelihood
  objective <- function(a) {
    # message("z:")
    # print(str(params$z / sum(params$z)))
    - sum( params$z * lpdf_transform(a) )
  }

  a0 <- mparam_rev_transform(params$theta);
  opt <- optim(a0, objective, method="L-BFGS-B", lower=-3, upper=3);
  theta <- mparam_transform(opt$par);

  lapply(seq_len(K), function(k) {
    list(mu = theta$mu[k], lambda = theta$lambda[k])
  })
}

# update_theta_beta <- function(C, v, params, hparams) {
#   Gamma <- params$Gamma;
#   theta <- params$theta;
#   K <- dim(Gamma)[3]
#   lapply(seq_len(K), function(k) {
#     weights <- colSums(C * Gamma[, , k, drop = FALSE][, , 1])
#     total <- sum(weights)
#     if (total <= 0) {
#       return(theta[[k]])
#     }

#     mu_hat <- sum(weights * v) / total
#     var_hat <- sum(weights * (v - mu_hat)^2) / total
#     lambda_hat <- if (var_hat <= 0) {
#       hparams$lambda_bounds[2]
#     } else {
#       mu_hat * (1 - mu_hat) / var_hat - 1
#     }
#     lambda_hat <- max(hparams$lambda_bounds[1], min(hparams$lambda_bounds[2], lambda_hat))
#     mu_hat <- max(hparams$mu_eps, min(1 - hparams$mu_eps, mu_hat))

#     list(mu = mu_hat, lambda = lambda_hat)
#   })
# }

set.seed(2)
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
target <- true_w %*% pmf_mat;
target <- target / rowSums(target);
C <- t(apply(target, 1, function(p) as.vector(rmultinom(1, counts_total, p))))

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

pmf_mat_hat <- t(matrix(unlist(fit$params$f), nrow=ncol(C), ncol=ncol(true_w)));
plot(pmf_mat, pmf_mat_hat)
cor(c(pmf_mat), c(pmf_mat_hat))

target.hat <- fit$params$w %*% pmf_mat_hat;
rowSums(target.hat)
target.hat
target

plot(target, target.hat)
cor(t(target), t(target.hat))
cor(c(target), c(target.hat))

i <- 3;
plot(v, target[i, ])
points(v, target.hat[i, ], col="blue")

cat("Mixture weights:\n")
print(round(fit$params$w, 3))
cat("Component parameters:\n")
print(lapply(fit$params$theta, function(th) list(mu = round(th$mu, 3), lambda = round(th$lambda, 3))))
