cohen.kappa <- function(xtab)
{

## xtab = table(p/a response data,glm.object$fitted>0.5)

# computes Cohen's kappa and variance estimates, 95% CI
# from a symmetric agreement table xtab
# return value is a list with elements kap and vark
# See Bishop et al. Disc. Mult. Analysis pp. 395-397 1975
# written by Lynd D. Bacon (USA)
	ci <- vector(length = 2)
	kap <- NULL
	vark <- NULL
	totn <- NULL
	if(nrow(xtab) != ncol(xtab)) {
		cat("\n freq. table not symmetric.\n")
		return(kap, vark, totn)
	}
# compute obs props:
	totn <- sum(xtab)
	obs <- xtab/totn	# calc marginals, expected diag props.
	pi <- apply(xtab, 1, sum)/totn
	pj <- apply(xtab, 2, sum)/totn
	exp <- outer(pi, pj)
	theta.2 <- sum(diag(exp))
	theta.1 <- sum(diag(obs))
	kap <- (theta.1 - theta.2)/(1 - theta.2)
	theta.3 <- diag(obs) %*% (pi + pj)
	theta.4 <- sum(obs * outer(pi, pj, FUN = "+")^2)
	# cat("\n thetas 1,2,3,4 \n", theta.1, theta.2, theta.3, theta.4, "\n")
# calc the var est
	vark <- (theta.1 * (1 - theta.1))/(1 - theta.2)^2
	vark <- vark + (2 * (1 - theta.1) * (2 * theta.1 * theta.2 - theta.3))/(1 - theta.2)^3
	vark <- vark + ((1 - theta.1)^2 * (theta.4 - 4 * theta.2^2))/(1 - theta.2)^4
	vark <- as.vector(vark/totn)
	ci[1] <- kap - 1.96 * sqrt(vark)
	ci[2] <- kap + 1.96 * sqrt(vark)
	return(list(kap=kap, vark=vark, totn=totn, ci=ci))
}
