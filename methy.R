library(io)
library(devtools)
load_all("pmfmix")

# data <- t(qread("data/ccoc-ts_methy_beta_cleaned_dist.rds"));
pheno0 <- qread("data/ccoc-ts_sample-info_stage3.tsv");
pheno0$cluster[pheno0$cluster == "mixture"] <- "hypomethylated";
pheno0$cluster[which(pheno0$sample_type == "normal")] <- "normal";

idx <- match(rownames(data[-1, ]), pheno0$sample_id);
pheno <- pheno0[idx, ];
all(pheno$sample_id == rownames(data[-1, ]), na.rm=TRUE)

# ---

beta_pmf <- function(v, theta) {
  dbeta(v, theta$mu * theta$lambda, (1 - theta$mu) * theta$lambda);
}

initialize_theta_beta <- function(C, v, K, hparams) {
  mu_grid <- seq(min(v), max(v), length.out = K + 2)[2:(K+1)];
  lambda <- rgamma(K, 10, 1);
  lapply(seq_len(K), function(k) {
    list(mu = mu_grid[k], lambda = lambda[k])
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
  lpdf_transform <- function(a) {
    theta <- mparam_transform(a);
    lp <- with(theta, unlist(lapply(v,
      # mixture of beta distributions
      function(x) log(t(params$W) * dbeta(x, mu*lambda, (1 - mu)*lambda))
    )));
    # W^T is K by N,  dbeta(x, ...) is K   ->  each item is K by N
    # output is K by N by J; need N by J by K
    lp <- aperm(array(lp, c(K, N, J)), c(2, 3, 1))
    lp
  }

  # negative log likelihood
  objective <- function(a) {
    - sum( params$Z * lpdf_transform(a) )
  }

  a0 <- mparam_rev_transform(params$theta);
  opt <- optim(a0, objective, method="L-BFGS-B", lower=-3, upper=3);
  theta <- mparam_transform(opt$par);

  lapply(seq_len(K), function(k) {
    list(mu = theta$mu[k], lambda = theta$lambda[k])
  })
}

# ---

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
    alpha = rep(1, K)
  ),
  control = list(nstart = 3, niter = 25, abstol = 1e-4),
  fixed = list(W = NULL, theta = NULL, Gamma = NULL),
  verbose = TRUE
)

predicted.pdfs <- t(vapply(fit$params$theta, function(th) beta_pmf(v, th), numeric(length(v))))
predicted.pdfs <- fit$params$W %*% predicted.pdfs
predicted.pdfs <- predicted.pdfs / rowSums(predicted.pdfs)
yhat2 <- predicted.pdfs;

sum( (yhat2 - target.pdfs)^2 )

fit
head(fit$params$W)
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

dim(fit$params$W)
rowSums(fit$params$W)

library(mmalign)
qdraw(
	pca_plot(t(fit$params$W), pheno=pheno, aes(colour=cluster)),
	width = 6,
	file = insert(out.fn, tag=c("params", "pca"), ext="pdf")
)

qwrite(fit$params, insert(out.fn, tag="params", ext="rds"));

