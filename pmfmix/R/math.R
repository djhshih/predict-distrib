logistic <- function(x) {
	1 / (1 + exp(-x))
}

logit <- function(x) {
  log(x) - log(1 - x)
}
