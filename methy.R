library(io)
source("R/pmfmix.R")

# data <- t(qread("data/ccoc-ts_methy_beta_cleaned_dist.rds"));
pheno0 <- qread("data/ccoc-ts_sample-info_stage3.tsv");
pheno0$cluster[pheno0$cluster == "mixture"] <- "hypomethylated";
pheno0$cluster[which(pheno0$sample_type == "normal")] <- "normal";

idx <- match(rownames(data[-1, ]), pheno0$sample_id);
pheno <- pheno0[idx, ];
all(pheno$sample_id == rownames(data[-1, ]), na.rm=TRUE)

# ---

beta_pmf <- function(v, theta) {
  probs <- dbeta(v, theta$mu * theta$lambda, (1 - theta$mu) * theta$lambda);
  probs / sum(probs)
}

initialize_theta_beta <- function(C, v, K, hparams) {
  mu_grid <- seq(min(v), max(v), length.out = K + 2)[2:(K+1)];
  lambda <- rgamma(K, 10, 1);
  lapply(seq_len(K), function(k) {
    list(mu = mu_grid[k], lambda = lambda[k])
  })
}

# FIXME these updates are based on heuristics
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
    lambda_hat <- if (var_hat <= 0) hparams$lambda_bounds[2] else mu_hat * (1 - mu_hat) / var_hat - 1
    mu_hat <- max(hparams$mu_eps, min(1 - hparams$mu_eps, mu_hat))
    lambda_hat <- max(hparams$lambda_bounds[1], min(hparams$lambda_bounds[2], lambda_hat))

    list(mu = mu_hat, lambda = lambda_hat)
  })
}

v <- data[1, ];
target.pdfs <- data[-1, ];
target.pdfs <- target.pdfs / rowSums(target.pdfs);
N <- 915452;
counts <- round(target.pdfs * N);

dim(counts)
dim(target.pdfs)

# ---

K <- 30
fit <- pmfmix(
  C = counts,
  v = v,
  K = K,
  f = beta_pmf,
  initialize_theta = initialize_theta_beta,
  update_theta = update_theta_beta,
  hparams = list(
    alpha = rep(1, K),
    lambda_bounds = c(2, 500),
    mu_eps = 1e-6
  ),
  control = list(nstart = 3, niter = 25, abstol = 1e-4),
  fixed = list(w = NULL, theta = NULL, Gamma = NULL),
  verbose = TRUE
)

predicted.pdfs <- t(vapply(fit$params$theta, function(th) beta_pmf(v, th), numeric(length(v))))
predicted.pdfs <- fit$params$w %*% predicted.pdfs
predicted.pdfs <- predicted.pdfs / rowSums(predicted.pdfs)
yhat2 <- predicted.pdfs;

sum( (yhat2 - target.pdfs)^2 )

fit
head(fit$params$w)
lapply(fit$params$theta, unlist)

out.fn <- filename("ccoc-ts", path="out", tag=c("methy", "beta", "pmfmix"));


mses <- rowSums((target.pdfs - yhat2)^2);
idx <- order(mses, decreasing=TRUE);

i <- idx[200];
plot(v, target.pdfs[i, ], type="l")
lines(v, counts[i, ] / rowSums(counts[i, , drop=FALSE]), col="red")
lines(v, yhat2[i, ], col="blue")

# @params xs  1 by M
# @params ys  N by M
distrib_plots <- function(xs, ys, main="pdf", type="l") {
	plot(xs, apply(ys, 2, max), main=main, xlab="x", ylab="f(x)", type=type)
	for (i in 1:nrow(ys)) {
		lines(xs, ys[i, ], main=main, col=i, type=type)
	}
}

qdraw(
	{
		par(mfrow=c(2, 1))
		distrib_plots(v, target.pdfs, main="observed pdf")
		distrib_plots(v, yhat2, main="predicted pdf")
	},
	height = 10,
	file = insert(out.fn, tag=c("pdf", "observed-vs-predicted"), ext="pdf")
)

sum( (yhat2 - target.pdfs)^2 )

dim(fit$params$w)
rowSums(fit$params$w)

library(mmalign)
qdraw(
	pca_plot(t(fit$params$w), pheno=pheno, aes(colour=cluster)),
	width = 6,
	file = insert(out.fn, tag=c("params", "pca"), ext="pdf")
)

qwrite(fit$params, insert(out.fn, tag="params", ext="rds"));

