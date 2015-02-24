SQRT2 = 1.41421356237309504880168872421

sgn = function(x) {
	x[x >= 0] = 1
	x[x < 0] = -1
	return(x)
}

dsgt = function(x, mu = 0, sigma = SQRT2, lambda = 0, p = 2, q = 100, mean.cent = TRUE, log = FALSE) {
	n = max(length(x), length(mu), length(sigma), length(lambda), length(p), length(q))
	x = rep(x, length.out=n)
	mu = rep(mu, length.out=n)
	sigma = rep(sigma, length.out=n)
	lambda = rep(lambda, length.out=n)
	p = rep(p, length.out=n)
	q = rep(q, length.out=n)
	if(mean.cent[1L]) x = x + (2*sigma*lambda*q^(1/p)*beta(2/p,q-1/p))/beta(1/p,q)
	if(!log[1L]) return(p/(2*sigma*q^(1/p)*beta(1/p,q)*(1+abs(x-mu)^p/(q*sigma^p*(1+lambda*sgn(x-mu))^p))^(q+1/p)))
	return(log(p)-log(2)-log(sigma)-log(q)/p-lbeta(1/p,q)-(1/p+q)*log(1+abs(x-mu)^p/(q*sigma^p*(1+lambda*sgn(x-mu))^p)))
}

psgt = function(quant, mu = 0, sigma = SQRT2, lambda = 0, p = 2, q = 100, mean.cent = TRUE, lower.tail = TRUE, log.p = FALSE) {
	n = max(length(quant), length(mu), length(sigma), length(lambda), length(p), length(q))
	quant = rep(quant, length.out=n)
	mu = rep(mu, length.out=n)
	sigma = rep(sigma, length.out=n)
	lambda = rep(lambda, length.out=n)
	p = rep(p, length.out=n)
	q = rep(q, length.out=n)
	if(mean.cent[1L]) quant = quant + (2*sigma*lambda*q^(1/p)*beta(2/p,q-1/p))/beta(1/p,q)
	quant = quant - mu
	flip = quant > 0
	lambda[flip] = -lambda[flip]
	quant[flip] = -quant[flip]
	out = (1-lambda)/2+(lambda-1)/2*pbeta(1/(1+q*(sigma*(1-lambda)/(-quant))^p),1/p,q)
	out[flip] = 1 - out[flip]
	if(!lower.tail[1L]) out = 1 - out
	if(log.p[1L]) out = log(out)
	return(out)
}

qsgt = function(prob, mu = 0, sigma = SQRT2, lambda = 0, p = 2, q = 100, mean.cent = TRUE, lower.tail = TRUE, log.p = FALSE) {
	n = max(length(prob), length(mu), length(sigma), length(lambda), length(p), length(q))
	prob = rep(prob, length.out=n)
	mu = rep(mu, length.out=n)
	sigma = rep(sigma, length.out=n)
	lambda = rep(lambda, length.out=n)
	p = rep(p, length.out=n)
	q = rep(q, length.out=n)
	if(log.p[1L]) prob = exp(prob)
	if(!lower.tail[1L])	prob = 1 - prob
	flip = prob > (1-lambda)/2
	prob[flip] = 1 - prob[flip]
	lam = lambda
	lam[flip] = -lam[flip]
	out = sigma*(lam-1)*(1/(q*qbeta(1-2*prob/(1-lam),1/p,q))-1/q)^(-1/p)
	out[flip] = -out[flip]
	out = out + mu
	if(mean.cent[1L]) out = out - (2*sigma*lambda*q^(1/p)*beta(2/p,q-1/p))/beta(1/p,q)
	return(out)
}

rsgt = function(n, mu = 0, sigma = SQRT2, lambda = 0, p = 2, q = 100, mean.cent = TRUE) {
	if(length(n) > 1) n = length(n)
	mu = rep(mu, length.out=n)
	sigma = rep(sigma, length.out=n)
	lambda = rep(lambda, length.out=n)
	p = rep(p, length.out=n)
	q = rep(q, length.out=n)
	qsgt(runif(n), mu, sigma, lambda, p, q, mean.cent)
}