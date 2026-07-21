# x is in real space
# y is in constrained space

logistic <- function(x) {
	1 / (1 + exp(-x))
}

logit <- function(y) {
	log(y) - log(1 - y)
}

softmax <- function(x) {
	exp(x) / sum(exp(x))
}

inv_softmax <- function(y) {
	log(y) - mean(log(y))
}
