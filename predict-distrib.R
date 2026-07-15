# Predict a bounded continuous distribution using a mixture of betas

softmax <- function(x) {
	exp(x) / sum(exp(x))
}

logistic <- function(x) {
	1 / (1 + exp(-x))
}

param_transform <- function(a) {
	K <- length(a) / 3;
	theta <- c(softmax(a[1:K]), logistic(a[(K+1):(2*K)]), exp(a[(2*K+1):length(a)]));
	w <- theta[1:K];
	mu <- theta[(K+1):(2*K)];
	lambda <- theta[(2*K+1):(3*K)];

	list(w = w, mu = mu, lambda = lambda)
}

# Transform activities a to a cumulative distribution fucntion that is evaluated at xs
pdf_transform <- function(a, xs) {
	params <- param_transform(a);

	ys <- with(params, unlist(lapply(xs,
		# mixture of beta distributions
		function(x) sum(w * dbeta(x, mu * lambda, (1 - mu)*lambda))
	)));
	ys / sum(ys)
}

J <- 30;
a <- rnorm(J, sd=0.5);
eps <- 1e-6;
xs <- seq(0+eps, 1-eps, length.out=100);

yhat <- pdf_transform(a, xs);

plot(xs, yhat, type="l", main="predicted cdf")
plot(xs[-1], diff(yhat), type="l", main="predicted pdf")

# ---

# Use numeric optimization to fit a target CDF using the
# pdf_transformation function

target.pdf <- 0.1 * dbeta(xs, 10, 2) + 0.2 * dbeta(xs, 3, 35);
target.pdf <- target.pdf / sum(target.pdf);

N <- 10000;
samp <- sample.int(length(xs), N, replace=TRUE, prob=target.pdf);
counts <- hist(xs[samp], breaks=c(0, xs), plot=FALSE)$counts;

plot(xs, target.pdf, type="l")
points(xs, counts/sum(counts))

# fit ground-truth pdf

objective_pdf <- function(a) {
	sum((target.pdf - pdf_transform(a, xs))^2)
}

a0 <- rnorm(J);
opt <- optim(a0, objective_pdf, method="L-BFGS-B");

yhat <- pdf_transform(opt$par, xs);

plot(xs, target.pdf, type="l")
lines(xs, yhat, col="royalblue3");

# fit count data

pseudo.count <- 1;

# negative unnormalized log multinomial distribution
objective_ulmultinom <- function(a) {
	- sum( (counts + pseudo.count - 1) * log(pdf_transform(a, xs)) )
}

opt2 <- optim(a0, objective_ulmultinom, method="L-BFGS-B");
yhat2 <- pdf_transform(opt2$par, xs);

param_transform(opt2$par)

plot(xs, target.pdf, type="l")
lines(xs, yhat2, col="royalblue3");

# ---

library(io)
data <- qread("data/ccoc-ts_methy_beta_cleaned_dist.rds");

# fit on counts of one real sample

xs <- data[, 1];
target.pdf <- data[, 2];
target.pdf <- target.pdf / sum(target.pdf);
N <- 915452;
counts <- round(target.pdf * N);

a0 <- rnorm(J);
opt2 <- optim(a0, objective_ulmultinom, method="L-BFGS-B", control=list(maxit=1e3));
opt2
yhat2 <- pdf_transform(opt2$par, xs);

param_transform(opt2$par)

plot(xs, target.pdf, type="l")
lines(xs, yhat2, col="royalblue3");

