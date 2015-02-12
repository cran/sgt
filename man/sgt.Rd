\name{sgt}
\alias{SGT}
\alias{sgt}
\alias{dsgt}
\alias{psgt}
\alias{qsgt}
\alias{rsgt}

\title{The Skewed Generalized T Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the skewed generalized t distribution.
}
\usage{
dsgt(x, mu = 0, sigma = SQRT2, lambda = 0, p = 2, q = 100, 
mean.cent = TRUE, log = FALSE)
psgt(quant, mu = 0, sigma = SQRT2, lambda = 0, p = 2, q = 100, 
mean.cent = TRUE, lower.tail = TRUE, log.p = FALSE)
qsgt(prob, mu = 0, sigma = SQRT2, lambda = 0, p = 2, q = 100, 
mean.cent = TRUE, lower.tail = TRUE, log.p = FALSE)
rsgt(n, mu = 0, sigma = SQRT2, lambda = 0, p = 2, q = 100, 
mean.cent = TRUE)
}
\arguments{
  \item{x, quant}{vector of quantiles.}
  \item{prob}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{mu}{vector of parameters. Note that if \code{mean.cent == TRUE}, \code{mu} is the mean of the distribution. Otherwise, \code{mu} is the mode of the distribution.}
  \item{sigma}{vector of variance parameters. The default is \code{sqrt(2)}. The variance of the distribution increases as \code{sigma} increases. Must be strictly positive.}
  \item{lambda}{vector of skewness parameters. Note that \code{-1 < lambda < 1}. If \code{lambda < 0}, the distribution is skewed to the left. If \code{lambda > 0}, the distribution is skewed to the right. If \code{lambda = 0}, then the distribution is symmetric.}
  \item{p, q}{vector of parameters. Smaller values of \code{p} and \code{q} result in larger values for the kurtosis of the distribution. Note that \code{p > 0}, \code{q > 0}, and \code{p*q > 1} must hold.}
  \item{mean.cent}{logical; if TRUE, \code{mu} is the mean of the distribution, otherwise \code{mu} is the mode of the distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.}
}
\value{
  \code{dsgt} gives the density,
  \code{psgt} gives the distribution function,
  \code{qsgt} gives the quantile function, and
  \code{rsgt} generates random deviates.
  
  The length of the result is determined by \code{n} for
  \code{rsgt}, and is the maximum of the lengths of the
  numerical arguments for the other functions.  
  
  The numerical arguments other than \code{n} are recycled to the
  length of the result.  Only the first elements of the logical
  arguments are used.
  
  \code{sigma <= 0}, \code{lambda <= -1}, \code{lambda >= 1}, \code{p <= 0}, \code{q <= 0}, and \code{p*q <= 1} are errors and return \code{NaN}.
}
\author{Carter Davis, \email{carterdavis@byu.edu}}
\details{
  If \code{mu}, \code{sigma}, \code{lambda}, \code{p}, or \code{q} are not specified they assume the default values of \code{mu = 0}, \code{sigma} is approximately \code{sqrt(2)}, \code{lambda = 0}, \code{p = 2}, and \code{q = 100}. These default values approximate a standard normal. As \code{q} approaches infinity with the other default values fixed, this distribution converges to the standard normal distribution.

  The skewed generalized t distribution has density
  \deqn{ 	
  	    f(x) = \frac{p}{2 \sigma q^{1/p} Beta( \frac{1}{p}, q) ( \frac{| x-\mu+ m | ^p}{q \sigma ^p (\lambda sign(x-\mu+ m) + 1)^p} + 1)^{\frac{1}{p}+q}}
  }
  where \eqn{m} is the mean of the distribution, \eqn{sign(x) = 1} for \eqn{x >= 0}, and \eqn{sign(x) = -1} for \eqn{x < 0}. If \code{mean.cent = TRUE}, then \eqn{m} is set to zero. The h\emph{th} moment of the distribution is defined for \eqn{p*q > h} to be 
  \deqn{
  	    \frac{\sigma^h q^{\frac{h}{p}} Beta(\frac{h+1}{p},q-\frac{h}{p})}{2 Beta(\frac{1}{p},q)} ((1+\lambda)^{h+1}+(-1)^h (1-\lambda)^{h+1})
  }
  
  Note that for large values of \code{sigma}, both the \code{psgt} and \code{qsgt} algortithms can be overwhelmed and result in an error.
}
\seealso{
  \link{Distributions} for other standard distributions which are special cases of the skewed generalized t distribution, including \code{\link{dt}} for the t distribution, \code{\link{dnorm}} for the normal distribution, and \code{\link{dunif}} for the uniform distribution. Other special cases of the skewed generalized t distribution include the generalized t distribution in the \code{gamlss.dist} package, the skewed t distribution in the \code{skewt} package, the exponential power distribution (also known as the generalized error distribution) in the \code{normalp} package, and the Laplace distribution in the \code{rmutil} package. Also see \code{\link{beta}} for the beta function.
}
\source{
  For \code{psgt}, based on

  the GNU Scientific Library (GSL) QAG and QAGIL quadrature algorithms. This provides accurate results for at least 10 digits.
  https://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration.html 

  For \code{qnorm}, based on

  the GNU Scientific Library (GSL) Steffensen root finding algorithm. The algorithm finds the root of 
  \deqn{
  	F(x) - p
  	}
  where \eqn{F(x)} is the distribution function of the skewed generalized t distribution, and \eqn{p} is the given probability. This algorithm provides precise results for at least 9 digits past the decimal.
  https://www.gnu.org/software/gsl/manual/html_node/One-dimensional-Root_002dFinding.html

  For \code{rnorm}, the algorithm employs a simple slice sampler to sample the skewed generalized t distribution.
}
\references{
	
  Besag, J., and Green, P. J. (1993)
  "Spatial Statistics and Bayesian Computation"
  \emph{Journal of the Royal Statistical Society}, volume 55, 25-37.
  
  Hansen, C., McDonald, J. B., and Newey, W. K. (2010)
  "Instrumental Variables Regression with Flexible Distributions"
  \emph{Journal of Business and Economic Statistics}, volume 28, 13-25.
  
  Kerman, S. C., and McDonald, J. B. (2012)
  "Skewness-Kurtosis Bounds for the Skewed Generalized T and Related Distributions"
  \emph{Statistics and Probability Letters}, volume 83, 2129-2134.
	
  Theodossiou, Panayiotis (1998)
  "Financial Data and the Skewed Generalized T Distribution"
  \emph{Management Science}, volume 44, 1650-1661.
  
}
\examples{
require(graphics)

### This shows how default values of the SGT allows
### it to approximate a standord normal distribution
x = seq(-4,4,by=0.05)
plot(x, dnorm(x), type='l')
lines(x, dsgt(x), col='blue')

### This shows how the p and q parameters change the distribution
plot(x, dsgt(x, p = 2, q = 2), type='l')
lines(x, dsgt(x, p = 5, q = 1000), col='blue')

psgt(0)
psgt(0, lambda = -0.5)
}
\keyword{distribution}