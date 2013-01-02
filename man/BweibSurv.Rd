\name{BweibSurv}
\alias{BweibSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to fit a Bayesian parametric (Weibull) regression model to univariate survival data.
}
\description{
The function to fit a Bayesian parametric (Weibull) regression model to univariate survival data.
}
\usage{
BweibSurv(survObj, priorPara, initial, num.reps, thin = 1, chain = 1,
			save = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{survObj}{
				The list containing observed data from \code{n} subjects; 
				\code{y}, \code{delta}, \code{x}
				}
\item{priorPara}{
				The list containing prior parameter values; 
					\code{a}, \code{b}, \code{c}, \code{d}
					}
					
\item{initial}{
				The list containing the starting values of the parameters;
				\code{beta.ini}, \code{alpha.ini}, \code{kappa.ini}
				}
\item{num.reps}{
the number of random scans of the chain
}
  \item{thin}{
thinning
}
  \item{chain}{
the numeric name of chain in the case when running multiple chains.
}
  \item{save}{
frequency of storing the results in .Rdata file. 
For example, by setting "save = 1000", the algorithm saves the results every 1000 iterations.
}

}
\details{
\tabular{ll}{
\code{y} \tab a vector of \code{n} times to the event \cr
\code{delta} \tab a vector of \code{n} censoring indicators for the event time (1=event occurred, 0=censored) \cr
\code{x} \tab covariate matrix, \code{n} observations by \code{p} variables\cr
\code{a} \tab the shape parameter of conjugate gamma prior for \eqn{\alpha}\cr
\code{b} \tab the rate parameter of conjugate gamma prior for \eqn{\alpha}\cr
\code{c} \tab the shape parameter of conjugate gamma prior for \eqn{\kappa}\cr
\code{d} \tab the rate parameter of conjugate gamma prior for \eqn{\kappa}\cr
\code{beta.ini} \tab the starting values for \eqn{\beta}\cr
\code{alpha.ini} \tab the starting values for \eqn{\alpha}\cr
\code{kappa.ini} \tab the starting values for \eqn{\kappa}\cr
}
}

\value{
\code{BweibSurv} returns an object of class \code{BweibSurv} \cr
\item{move.ind}{selected moves in the MCMC sampling}
\item{beta.p}{posterior samples for \eqn{\beta}}
\item{alpha.p}{posterior samples for \eqn{\alpha}}
\item{kappa.p}{posterior samples for \eqn{\kappa}}
\item{accept.beta}{the number of acceptance in sampling \eqn{\beta}}
\item{accept.alpha}{the number of acceptance in sampling \eqn{\alpha}}
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2013). 
Bayesian Semi-parametric Analysis of Semi-competing Risks Data: 
Estimating Readmission Rates among Pancreatic Cancer Patients, available from the authors upon request.
}
\author{
Kyu Ha Lee, Sebastein Haneuse, Deborah Schrag, and Francesca Dominici \cr
}
\note{
If the prespecified value of \code{save} is less than that of \code{num.reps}, the results are saved
as \code{.Rdata} file under the directory \code{working directory/mcmcOutcome}. 
Trace plots of posterior samples for the parameters are also provided, with the first half taken as burn-in.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

	### simulating survival data ###
	
	set.seed(123456)
	
	n = 2000
	p = 2

	x	= matrix(0, n, p)	
	x[,1]	= rnorm(n, 0, 2)	
	x[,2]	= sample(c(0, 1), n, replace = TRUE)
			
	
	# - true value of parameter
	
	beta.true = c(0.5, 0.5)
	
	alpha.true = 1.5		
	kappa.true = 0.02

	cens <- c(30, 40)		

	simData <- SimSurv(x, beta.true, alpha.true, kappa.true, cens) 

   				  
	# Set prior parameters and initial values
	
	survObj	<- priorPara	<- initial	<- list()
		
	y <- survObj$y <- simData$y		
	delta <- survObj$delta <- simData$delta
	survObj$x	<- x

	priorPara$a			<- 0.5	# prior parameters for alpha
	priorPara$b			<- 0.01
	
	priorPara$c			<- 0.5	# prior parameters for kappa
	priorPara$d			<- 0.05
		
	initial$beta.ini		<- rep(1, p)	
	initial$alpha.ini		<- runif(1, 0, 2)
	initial$kappa.ini		<- 1
	
	# fitting parametric regression model to semi-competing risks data
	# set 'num.reps' to a larger number such as 500000
	
	fit.BweibSurv <- BweibSurv(survObj, priorPara, initial, num.reps = 100,
	 					thin = 10, chain = 1, save = 100)
	

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ parametric analysis }
\keyword{ univariate survival analysis }% __ONLY ONE__ keyword per line
