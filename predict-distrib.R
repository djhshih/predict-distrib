# Predict a bounded continuous distribution using a mixture of betas

softmax <- function(x) {
	exp(x) / sum(exp(x))
}

# Transform activities a to a cumulative distribution fucntion that is evaluated at xs
cdf_transform <- function(a, xs) {
	K <- length(a) / 3;
	theta <- c(softmax(a[1:K]), exp(a[(K+1):J]));
	pi <- theta[1:K];
	alpha <- theta[(K+1):(2*K)];
	beta <- theta[(2*K+1):(3*K)];

	unlist(lapply(xs,
		# mixture of beta distributions
		function(x) sum(pi * pbeta(x, alpha, beta))
	))
}

J <- 30;
a <- rnorm(J, sd=0.5);
xs <- seq(0, 1, by=0.01);

yhat <- cdf_transform(a, xs);

plot(xs, yhat, type="l", main="predicted cdf")
plot(xs[-1], diff(yhat), type="l", main="predicted pdf")

# ---

# Use numeric optimization to fit a target CDF using the
# cdf_transformation function

target.pdf <- 0.1 * dnorm(xs, 0.1, 0.2) + 0.2 * dnorm(xs, 0.95, 0.1);
target.cdf <- cumsum(target.pdf) / sum(target.pdf);

plot(xs, target.pdf, type="l")

objective <- function(a) {
	sum((target.cdf - cdf_transform(a, xs))^2)
}

a0 <- rnorm(J);
opt <- optim(a0, objective, method="CG");

yhat <- cdf_transform(opt$par, xs);

plot(xs, target.cdf, type="l")
lines(xs, yhat, col="royalblue3");

