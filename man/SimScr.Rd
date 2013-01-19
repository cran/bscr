\name{SimScr}
\alias{SimScr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to simulate semi-competing risks data under Weibull model
}
\description{
The function to simulate semi-competing risks data under Weibull model
}
\usage{
SimScr(x, beta1.true, beta2.true, beta3.true, alpha1.true, alpha2.true, 
		alpha3.true, kappa1.true, kappa2.true, kappa3.true, gamma.true, cens)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Covariate matrix, \code{n} observations by \code{p} variables
}
  \item{beta1.true}{
True value for \eqn{\beta_1}
}
  \item{beta2.true}{
True value for \eqn{\beta_2}
}
  \item{beta3.true}{
True value for \eqn{\beta_3}
}
  \item{alpha1.true}{
True value for \eqn{\alpha_1}
}
  \item{alpha2.true}{
True value for \eqn{\alpha_2}
}
  \item{alpha3.true}{
True value for \eqn{\alpha_3}
}
  \item{kappa1.true}{
True value for \eqn{\kappa_1}
}
  \item{kappa2.true}{
True value for \eqn{\kappa_2}
}
  \item{kappa3.true}{
True value for \eqn{\kappa_3}
}
  \item{gamma.true}{
True value for \eqn{\gamma}
}

  \item{cens}{
A vector with two numeric elements. The right censoring times are generated from Uniform(\eqn{cens[1]}, \eqn{cens[2]})
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{SimScr} returns an object of class \code{SimScr} \cr
\item{y1}{a vector of \code{n} times to the non-terminal event}
\item{y2}{a vector of \code{n} times to the terminal event}
\item{delta1}{a vector of \code{n} censoring indicators for the non-terminal event time (1=event occurred, 0=censored)}
\item{delta2}{a vector of \code{n} censoring indicators for the terminal event time (1=event occurred, 0=censored)}
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2013). 
Bayesian Semi-parametric Analysis of Semi-competing Risks Data: 
Estimating Readmission Rates among Pancreatic Cancer Patients, available from the authors upon request.
}
\author{
Kyu Ha Lee, Sebastien Haneuse, Deborah Schrag, and Francesca Dominici \cr
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
	### simulating semi-competing risks data ###
	
	set.seed(123456)
	
	n = 2000
	p = 2

	x	= matrix(0, n, p)	
	x[,1]	= rnorm(n, 0, 2)	
	x[,2]	= sample(c(0, 1), n, replace = TRUE)
			
				
	# - true value of parameter
	
	beta1.true = c(0.5, 0.5)
	beta2.true = c(0.5, 0.5)
	beta3.true = c(1, 1)
	
	theta.true = 2
	gamma.true = rgamma(n, shape = 1/theta.true, rate = 1/theta.true)

	alpha1.true = 1.5
	alpha2.true = 0.8
	alpha3.true = 0.8
		
	kappa1.true = 0.05
	kappa2.true = 0.05
	kappa3.true = 0.10
	
	cens <- c(30, 40)		

	simData <- SimScr(x, beta1.true, beta2.true, beta3.true, alpha1.true, 
    			alpha2.true, alpha3.true, kappa1.true, kappa2.true, kappa3.true, 
   				 gamma.true, cens) 
   	summary(simData)			 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line
