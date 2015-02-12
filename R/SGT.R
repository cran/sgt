SQRT2 = 1.41421356237309504880168872421

dsgt = function(x, mu = 0, sigma = SQRT2, lambda = 0, p = 2, q = 100, mean.cent = TRUE, log = FALSE) {
	outlength = max(c(length(x),length(mu),length(sigma),length(lambda),length(p),length(q)))
	if (!is.finite(mean.cent)) mean.cent = TRUE
	if (!is.finite(log)) log = TRUE
	fin = rep(TRUE, outlength)
	finList = list(x = x, mu = mu, sigma = sigma, lambda = lambda, p = p, q = q)
	for (i in 1:length(finList)) {
		fin = fin & rep(is.finite(finList[[i]]), outlength)
		eval(parse(text=paste(names(finList)[i], "[!is.finite(", names(finList)[i], ")] = 0")))
	}
	out = rep(0, outlength)
	retval = .C("dsgt",outval=as.double(out),as.integer(outlength),as.integer(fin),as.double(x),as.integer(length(x)),as.double(mu),as.integer(length(mu)),as.double(sigma),as.integer(length(sigma)),as.double(lambda),as.integer(length(lambda)),as.double(p),as.integer(length(p)),as.double(q),as.integer(length(q)),as.integer(mean.cent),as.integer(log))
	return(retval$outval)
}

psgt = function(quant, mu = 0, sigma = SQRT2, lambda = 0, p = 2, q = 100, mean.cent = TRUE, lower.tail = TRUE, log.p = FALSE) {
	outlength = max(c(length(quant),length(mu),length(sigma),length(lambda),length(p),length(q)))
	if (!is.finite(mean.cent)) mean.cent = TRUE
	if (!is.finite(lower.tail)) lower.tail = TRUE
	if (!is.finite(log.p)) log.p = TRUE
	fin = rep(TRUE, outlength)
	finList = list(quant = quant, mu = mu, sigma = sigma, lambda = lambda, p = p, q = q)
	for (i in 1:length(finList)) {
		fin = fin & rep(is.finite(finList[[i]]), outlength)
		eval(parse(text=paste(names(finList)[i], "[!is.finite(", names(finList)[i], ")] = 0")))
	}
	out = rep(0,outlength)
	retval = .C("psgt",outval=as.double(out),as.integer(outlength),as.integer(fin),as.double(quant),as.integer(length(quant)),as.double(mu),as.integer(length(mu)),as.double(sigma),as.integer(length(sigma)),as.double(lambda),as.integer(length(lambda)),as.double(p),as.integer(length(p)),as.double(q),as.integer(length(q)),as.integer(mean.cent),as.integer(lower.tail),as.integer(log.p))
	return(retval$outval)
}

qsgt = function(prob, mu = 0, sigma = SQRT2, lambda = 0, p = 2, q = 100, mean.cent = TRUE, lower.tail = TRUE, log.p = FALSE) {
	outlength = max(c(length(prob),length(mu),length(sigma),length(lambda),length(p),length(q)))
	out = rep(0,outlength)
	if (!is.finite(mean.cent)) mean.cent = TRUE
	if (!is.finite(lower.tail)) lower.tail = TRUE
	if (!is.finite(log.p)) log.p = TRUE
	fin = rep(TRUE, outlength)
	finList = list(prob = prob, mu = mu, sigma = sigma, lambda = lambda, p = p, q = q)
	for (i in 1:length(finList)) {
		fin = fin & rep(is.finite(finList[[i]]), outlength)
		eval(parse(text=paste(names(finList)[i], "[!is.finite(", names(finList)[i], ")] = 0")))
	}
	retval = .C("qsgt",outval=as.double(out),as.integer(outlength),as.integer(fin),as.double(prob),as.integer(length(prob)),as.double(mu),as.integer(length(mu)),as.double(sigma),as.integer(length(sigma)),as.double(lambda),as.integer(length(lambda)),as.double(p),as.integer(length(p)),as.double(q),as.integer(length(q)),as.integer(mean.cent),as.integer(lower.tail),as.integer(log.p))
	return(retval$outval)
}

rsgt = function(n, mu = 0, sigma = SQRT2, lambda = 0, p = 2, q = 100, mean.cent = TRUE) {
	if (!is.finite(mean.cent)) mean.cent = TRUE
    if (length(n) > 1) n = length(n)
	if (as.integer(n) <= 0) stop("'n' must be a postive integer")
	fin = rep(TRUE, n)
	finList = list(mu = mu, sigma = sigma, lambda = lambda, p = p, q = q)
	for (i in 1:length(finList)) {
		fin = fin & rep(is.finite(finList[[i]]), n)
		eval(parse(text=paste(names(finList)[i], "[!is.finite(", names(finList)[i], ")] = 0")))
	}
	burnin = 3
	out = rep(0,as.integer(n))
	unif = runif(2*burnin*as.integer(n))
	retval = .C("rsgt",outval=as.double(out),as.integer(n),as.integer(fin),as.integer(burnin),as.double(unif),as.double(mu),as.integer(length(mu)),as.double(sigma),as.integer(length(sigma)),as.double(lambda),as.integer(length(lambda)),as.double(p),as.integer(length(p)),as.double(q),as.integer(length(q)),as.integer(mean.cent))
	return(retval$outval)
}


