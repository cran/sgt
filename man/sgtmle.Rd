\name{sgtmle}
\alias{sgt.mle}
\alias{SGT.MLE}
\alias{print.MLE}
\alias{print.mult.MLE}

\title{Maximum Likelihood Estimation with the Skewed Generalized T Distribution}
\description{
  This function allows data to be fit to the skewed generalized t distribution using maximum likelihood estimation. This function uses the \code{maxLik} package to perform its estimations.
}
\usage{
sgt.mle(X.f, mu.f = mu ~ mu, sigma.f = sigma ~ sigma, 
lambda.f = lambda ~ lambda, p.f = p ~ p, q.f = q ~ q, 
data = parent.frame(), start, subset, method = 'BFGS', 
constraints = NULL, follow.on = FALSE, iterlim = 5000, \dots)
}
\arguments{
  \item{X.f}{A formula specifying the data, or the function of the data with parameters, that should be used in the maximisation procedure. \code{X} should be on the left-hand side and the right-hand side should be the data or function of the data that should be used.}
  \item{mu.f, sigma.f, lambda.f, p.f, q.f}{formulas including variables and
    parameters that specify the functional form of the parameters in the skewed generalized t log-likelihood function. \code{mu}, \code{sigma}, \code{lambda}, \code{p}, and \code{q} should be on the left-hand side of these formulas respectively.}
  \item{data}{an optional data frame in which to evaluate the variables in \code{formula} and \code{weights}.  Can also be a list or an environment.}
  \item{start}{a named list or named numeric vector of starting estimates for every parameter.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{method}{A list of the methods to be used. May include
   "NR" (for Newton-Raphson),
   "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno), 
   "BHHH" (for Berndt-Hall-Hall-Hausman), 
   "SANN" (for Simulated ANNealing), 
   "CG" (for Conjugate Gradients), 
   or "NM" (for Nelder-Mead).  
   Lower-case letters (such as "nr" for Newton-Raphson) are allowed.
   The default method is the "BFGS" method.}
  \item{constraints}{either \code{NULL} for unconstrained optimization
    or a list with two components.  The components may be either
    \code{eqA} and \code{eqB} for equality-constrained optimization
    \eqn{A \theta + B = 0}{A \%*\% theta + B = 0}; or \code{ineqA} and
    \code{ineqB} for inequality constraints \eqn{A \theta + B > 0}{A
      \%*\% theta + B > 0}.  More
       than one
       row in \code{ineqA} and \code{ineqB} corresponds to more than
       one linear constraint, in that case all these must be zero
       (equality) or positive (inequality constraints).}
  \item{follow.on}{logical; if TRUE, and there are multiple methods, then the last set of parameters from one method is used as the starting set for the next.}
  \item{iterlim}{If provided as a vector of the same length as \code{method}, gives the maximum number of iterations or function values for the corresponding method. If a single number is provided, this will be used for all methods.}
  \item{\dots}{further arguments that are passed to the selected maximisation routine in the \code{maxLik} package. See below for a non-exhaustive list of some further arguments that can be used.}
  }
\details{
  The parameter names are taken from \code{start}. If there is a name of a parameter or some data found on the right-hand side of one of the formulas but not found in \code{data} and not found in \code{start}, then an error is given.
  
  Below is a non-exhaustive list of further arguments that may be passed in to the \code{sgt.mle} function (see \code{maxLik} documentation for more details):
  \describe{
      \item{\code{fixed}}{
      
      parameters that should be fixed at their starting values:
      a vector of character strings indicating the names of the fixed
      parameters (parameter names are taken from argument \code{start}). May not be used in \code{BHHH} algorithm.}
  \item{\code{print.level}}{a larger number prints more working information.}
  \item{\code{tol}, \code{reltol}}{the absolute and relative convergence tolerance (see
    \code{\link{optim}}). May not be used in \code{BHHH} algorithm.}
  \item{\code{finalHessian}}{how (and if) to calculate the final Hessian.  Either
  \code{FALSE} (not calculate), \code{TRUE} (use analytic/numeric
  Hessian) or \code{"bhhh"}/\code{"BHHH"} for information equality
  approach.}
  \item{\code{parscale}}{A vector of scaling values for the parameters.
          Optimization is performed on 'par/parscale' and these should
          be comparable in the sense that a unit change in any element
          produces about a unit change in the scaled value. (see
    \code{\link{optim}}). May not be used in \code{BHHH} algorithm.}
    }
  Note that not all arguments may be used for every maximisation algorithm at this time. If multiple methods are supplied (i.e. \code{length(method) > 1}), all arguments are employed for each method (except for \code{iterlim}, which is allowed to vary for different methods).
  
If multiple methods are supplied, and some methods fail to initialise properly, a warning will be given. If every method fails to initialise, an error is given.
}
\value{
  If only one method is specified, \code{sgt.mle} returns a list of class \code{"MLE"}.
  If multiple methods are given, \code{sgt.mle} returns a list of class \code{"mult.MLE"} with each component containing the results of each maximisation procedure. Each component is a list of class \code{"MLE"}.
  A list of class \code{"MLE"} has the following components:
  \item{parameters}{the names of the given parameters taken from \code{start}}
  \item{maximum}{\code{fn} value at maximum (the last calculated value
    if not converged).}
  \item{estimate}{estimated parameter value.}
  \item{gradient}{vector, last gradient value which was calculated.  Should be
    close to 0 if normal convergence.}
  \item{gradientObs}{matrix of gradients at parameter value \code{estimate}
      evaluated at each observation (only if \code{grad} returns a matrix
      or \code{grad} is not specified and \code{fn} returns a vector).}
  \item{hessian}{Hessian at the maximum (the last calculated value if
    not converged).}
  \item{code}{return code:
    \itemize{
    \item{1}{ gradient close to zero (normal convergence).}
    \item{2}{ successive function values within tolerance limit (normal
      convergence).}
    \item{3}{ last step could not find higher value (probably not
      converged).  This is related to line search step getting too
      small, usually because hitting the boundary of the parameter
      space.  It may also be related to attempts to move to a wrong
      direction because of numerical errors.  In some cases it can be
      helped by changing \code{steptol}.}
    \item{4}{ iteration limit exceeded.}
    \item{5}{ Infinite value.}
    \item{6}{ Infinite gradient.}
    \item{7}{ Infinite Hessian.}
    \item{8}{ Successive function values withing relative tolerance
      limit (normal convergence).}
    \item{9}{ (BFGS) Hessian approximation cannot be improved because of
      gradient did not change.  May be related to numerical
      approximation problems or wrong analytic gradient.}
    \item{100}{ Initial value out of range.}
    }
  }
  \item{message}{ a short message, describing \code{code}.}
  \item{last.step}{ list describing the last unsuccessful step if
    \code{code=3} with following components:
    \itemize{
    \item{theta0}{ previous parameter value}
    \item{f0}{ \code{fn} value at \code{theta0}}
    \item{climb}{ the movement vector to the maximum of the quadratic approximation}
    }
  }
  \item{fixed}{logical vector, which parameters are constants.}
  \item{iterations}{number of iterations.}
  \item{type}{character string, type of maximization.}
  \item{constraints}{A list, describing the constrained optimization
    (\code{NULL} if unconstrained).  Includes the following components:
    \itemize{
      \item{type}{ type of constrained optimization}
      \item{outer.iterations}{ number of iterations in the constraints step}
      \item{barrier.value}{ value of the barrier function}
    }
  }
}
\author{Carter Davis, \email{carterdavis@byu.edu}}
\references{
  Henningsen, Arne and Toomet, Ott (2011). 
  "maxLik: A package for maximum likelihood estimation in R" 
  \emph{Computational Statistics} 26(3), 443-458. DOI 10.1007/s00180-010-0217-1.
}

\seealso{The \code{maxLik} package and its documentation. The \code{sgt.mle} simply uses its functions to maximize the skewed generalized t log-likelihood.
}
\examples{
### Showing how to fit a simple vector of data to the skewed 
### generalized t distribution. 
require(graphics)
require(stats)
set.seed(123456)
x = rt(100, df=10)
X.f = X ~ x
start = list(mu = 0, sigma = 2, lambda = 0, p = 2, q = 12)
result = sgt.mle(X.f = X.f, start = start, finalHessian = "BHHH")
sumResult = summary(result)
print(result)
coef(result)
print(sumResult)
### Note that the t distribution is a special case of the 
### skewed generalized t distribution

### Looking at the fit graphically
xvals = seq(-3,3,by=0.01)
plot(xvals, dt(xvals, df=10), type="l", col="blue")
lines(density(x))
mu = result$estimate[1]
sigma = result$estimate[2]
lambda = result$estimate[3]
p = result$estimate[4]
q = result$estimate[5]
lines(xvals, dsgt(xvals, mu = mu, sigma = sigma, lambda = lambda, p = p, q = q), col="red")
}
\keyword{optimize}
