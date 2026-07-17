# Predict bounded continuous distributions using a mixture of betas
# with global mean and scale parameters

softmax <- function(x) {
	exp(x) / sum(exp(x))
}

logistic <- function(x) {
	1 / (1 + exp(-x))
}

mparam_transform <- function(a, N, K) {
	w_end <- N * K;
	mu_start <- w_end + 1;
	mu_end <- w_end + K;
	lambda_start <- mu_end + 1;
	lambda_end <- mu_end + K;
	
	# w is N by K
	w <- t(apply(matrix(a[1:w_end], nrow=N), 1, softmax));
	mu <- logistic(a[mu_start:mu_end]);
	lambda <- exp(a[lambda_start:lambda_end]);

	list(w = w, mu = mu, lambda = lambda)
}

# Transform activities a to a probability mass function that is evaluated at xs
# return N x J matrix, where each row is a probability mass function
mpdf_transform <- function(a, xs, N, K) {
	params <- mparam_transform(a, N, K);

	ys <- with(params, matrix(unlist(lapply(xs,
		# mixture of beta distributions
		function(x) colSums(t(w) * dbeta(x, mu*lambda, (1 - mu)*lambda))
	)), nrow=N));
	# ys is N x J, where J = length(xs)

	ys / rowSums(ys)
}

N <- 10;
K <- 3;

J <- N*K + 2*K;
a <- rnorm(J, sd=0.5);
eps <- 1e-6;
M <- 100;
xs <- seq(0+eps, 1-eps, length.out=M);

params <- mparam_transform(a, N, K);
stopifnot(abs(rowSums(params$w) - 1) < 1e-9)

yhat <- mpdf_transform(a, xs, N, K);

plot(xs, yhat[1, ], type="l", main="predicted cdf")
for (i in 1:N) {
	lines(xs, yhat[i, ], type="l", main="predicted cdf")
}

# ---

# Use numeric optimization to fit a target CDF using the
# pdf_transformation function

target.pdfs <- matrix(0, nrow=N, ncol=M);
w <- runif(N);
w <- matrix(c(w, 1 - w), nrow=N);
mu <- c(10/12, 3/38);
lambda <- c(12, 38);
for (i in 1:N) {
	target.pdfs[i, ] <- unlist(lapply(xs, function(x) sum(w[i, ] * dbeta(x, mu*lambda, (1-mu)*lambda))));
}
target.pdfs <- target.pdfs / rowSums(target.pdfs);
rowSums(target.pdfs)

# @params xs  1 by M
# @params ys  N by M
distrib_plots <- function(xs, ys, main="pdf", type="l") {
	plot(xs, apply(ys, 2, max), main=main, xlab="x", ylab="f(x)", type=type)
	for (i in 1:N) {
		lines(xs, ys[i, ], main=main, col=i, type=type)
	}
}

distrib_plots(xs, target.pdfs)

# sample from target pdf

total.count <- 1000;
counts <- matrix(0, nrow=N, ncol=M);
for (i in 1:N) {
	samp <- sample.int(length(xs), total.count, replace=TRUE, prob=target.pdfs[i, ]);
	counts[i, ] <- hist(xs[samp], breaks=c(0, xs), plot=FALSE)$counts;
}

distrib_plots(xs, counts / rowSums(counts), type="b")


# fit ground-truth pdf


objective_pdf <- function(a) {
	sum((target.pdfs - mpdf_transform(a, xs, N, K))^2)
}

a0 <- rnorm(J);
opt <- optim(a0, objective_pdf, method="L-BFGS-B");
opt

yhat <- mpdf_transform(opt$par, xs, N, K);
stopifnot(dim(target.pdfs) == dim(yhat))

par(mfrow=c(2, 1))
distrib_plots(xs, target.pdfs, main="observed pdf")
distrib_plots(xs, yhat, main="predicted pdf")

# fit count data

# negative unnormalized log multinomial distribution
objective_ulmultinom <- function(a) {
	- sum( counts * log(mpdf_transform(a, xs, N, K)) )
}

opt2 <- optim(a0, objective_ulmultinom, method="L-BFGS-B");
opt2

yhat2 <- mpdf_transform(opt2$par, xs, N, K);
stopifnot(dim(counts) == dim(yhat2))

mparam_transform(opt2$par, N, K)

par(mfrow=c(2, 1))
distrib_plots(xs, target.pdfs, main="observed pdf")
distrib_plots(xs, yhat2, main="predicted pdf")

mean( (yhat - target.pdfs)^2 )
mean( (yhat2 - target.pdfs)^2 )

# ---

library(io)
data <- t(qread("data/ccoc-ts_methy_beta_cleaned_dist.rds"));
pheno0 <- qread("data/ccoc-ts_sample-info_stage3.tsv");
pheno0$cluster[pheno0$cluster == "mixture"] <- "hypomethylated";
pheno0$cluster[which(pheno0$sample_type == "normal")] <- "normal";

idx <- match(rownames(data[-1, ]), pheno0$sample_id);
pheno <- pheno0[idx, ];
all(pheno$sample_id == rownames(data[-1, ]), na.rm=TRUE)

# fit on counts of one real sample


xs <- data[1, ];
target.pdfs <- data[-1, ];
target.pdfs <- target.pdfs / rowSums(target.pdfs);
N <- 915452;
counts <- round(target.pdfs * N);

dim(counts)
dim(target.pdfs)

N <- nrow(target.pdfs);
K <- 10;
J <- N*K + 2*K;

a0 <- rnorm(J);
yhat0 <- mpdf_transform(a0, xs, N, K);

# NB: this takes several hours
opt2 <- optim(a0, objective_ulmultinom, method="L-BFGS-B", control=list(maxit=1e3));
opt2

# Do another round of optimization
a0 <- opt2$par;
opt2 <- optim(a0, objective_ulmultinom, method="L-BFGS-B", control=list(maxit=1e3));
opt2

yhat2 <- mpdf_transform(opt2$par, xs, N, K);

params <- mparam_transform(opt2$par, N, K);

out.fn <- filename("ccoc-ts", path="out", tag=c("methy", "beta"));

i <- 5;
plot(xs, target.pdfs[i, ], type="l")
lines(xs, counts[i, ] / rowSums(counts[i, , drop=FALSE]), col="red")
lines(xs, yhat2[i, ], col="blue")

qdraw(
	{
		par(mfrow=c(2, 1))
		distrib_plots(xs, target.pdfs, main="observed pdf")
		distrib_plots(xs, yhat2, main="predicted pdf")
	},
	height = 10,
	file = insert(out.fn, tag=c("pdf", "observed-vs-predicted"), ext="pdf")
)

mean( (yhat2 - target.pdfs)^2 )

library(mmalign)
qdraw(
	pca_plot(t(params$w), pheno=pheno, aes(colour=cluster)),
	width = 6,
	file = insert(out.fn, tag=c("mix-distrib-params", "pca"), ext="pdf")
)

qwrite(params, insert(out.fn, tag="mix-distrib-params", ext="rds"));

