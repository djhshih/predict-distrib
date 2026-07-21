library(io)
library(devtools)
load_all("pmfmix")

data <- t(qread("data/ccoc-ts_methy_beta_cleaned_dist.rds"));
pheno0 <- qread("data/ccoc-ts_sample-info_stage3.tsv");
pheno0$cluster[pheno0$cluster == "mixture"] <- "hypomethylated";
pheno0$cluster[which(pheno0$sample_type == "normal")] <- "normal";

rownames(data) <- sub("-02", "", rownames(data));
rownames(data) <- sub("R$", "", rownames(data));

idx <- match(rownames(data[-1, ]), pheno0$sample_id);
pheno <- pheno0[idx, ];
stopifnot(pheno$sample_id == rownames(data[-1, ]))

# ---

beta_lpmf <- function(v, theta) {
	with(theta, dbeta(v, mu * lambda, (1 - mu) * lambda, log=TRUE))
}

initialize_theta_beta <- function(C, v, K, hparams) {
	mu <- seq(min(v), max(v), length.out = K + 2)[2:(K+1)];
	lambda <- rgamma(K, 10, 1);
	list(mu = mu, lambda = lambda)
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
		c(logit(theta$mu), log(theta$lambda))
	}

	# Transform activities a to a probability mass function that is evaluated at xs
	# return N x M matrix, where each row is a probability mass function
	lwf_transform <- function(a) {
		theta <- mparam_transform(a);
		lp <- unlist(lapply(v,
			# mixture of beta distributions
			function(x) log(t(params$W)) + beta_lpmf(x, theta)
		));
		# W^T is K by N,  beta_lpmf(x, theta) is K   ->  each item is K by N
		# output is K by N by J; need N by J by K
		aperm(array(lp, c(K, N, J)), c(2, 3, 1))
	}

	# - E_Z[ log p(C, Z, theta) ]
	objective_q <- function(a) {
		- sum( params$Z * lwf_transform(a) );
	}

	a0 <- mparam_rev_transform(params$theta);
	opt <- optim(a0, objective_q, method="L-BFGS-B", lower=-10, upper=10);
	theta <- mparam_transform(opt$par);

	theta
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

set.seed(1234)

K <- 30
fit <- pmfmix(
	C = counts,
	v = v,
	K = K,
	lf = beta_lpmf,
	initialize_theta = initialize_theta_beta,
	update_theta = update_theta_beta,
	hparams = list(alpha = rep(1, K)),
	control = list(nstart = 3, niter = 25, abstol = 1e-4),
	verbose = TRUE
)

predicted.pdfs <- fit$params$W %*% exp(fit$params$lF);
rowSums(predicted.pdfs)

mean( (predicted.pdfs - target.pdfs)^2 )

fit
head(fit$params$W)

out.fn <- filename("ccoc-ts", path="out", tag=c("methy", "beta", "pmfmix"));


mses <- rowMeans((predicted.pdfs - target.pdfs)^2);
idx <- order(mses, decreasing=TRUE);

i <- idx[2];
plot(v, target.pdfs[i, ], type="l")
lines(v, counts[i, ] / rowSums(counts[i, , drop=FALSE]), col="red")
lines(v, predicted.pdfs[i, ], col="blue")

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
		distrib_plots(v, predicted.pdfs, main="predicted pdf")
	},
	height = 10,
	file = insert(out.fn, tag=c("pdf", "observed-vs-predicted"), ext="pdf")
)

mean( (predicted.pdfs - target.pdfs)^2 )

library(mmalign)
qdraw(
	pca_plot(t(fit$params$W), pheno=pheno, aes(colour=cluster)),
	width = 6,
	file = insert(out.fn, tag=c("params", "pca"), ext="pdf")
)

wsums <- colSums(fit$params$W);
idx <- order(wsums)[1:10];
qdraw(
	pca_plot(t(fit$params$W[, idx]), pheno=pheno, aes(colour=cluster)),
	width = 6,
	file = insert(out.fn, tag=c("top-w", "pca"), ext="pdf")
)

qwrite(fit$params, insert(out.fn, tag="params", ext="rds"));

