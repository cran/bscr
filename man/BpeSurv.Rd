\name{BpeSurv}
\alias{BpeSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to fit a Bayesian semi-parametric regression model to univariate survival data.
}
\description{
The function to fit a Bayesian semi-parametric regression model to univariate survival data.
}
\usage{
BpeSurv(survObj, priorPara, initial, num.reps, thin = 1, chain = 1,
		 save = 1000, RJ = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{survObj}{
				The list containing observed data from \code{n} subjects; 
				\code{y}, \code{delta}, \code{x}
				}
\item{priorPara}{
				The list containing prior parameter values; 
					\code{c.lam}, \code{C}, \code{a}, \code{b}, \code{alpha}, \code{delPert},
					\code{smax}, \code{s.propBI}
					}
					
\item{initial}{
				The list containing the starting values of the parameters;
				\code{beta.ini}
				\code{s}, \code{lambda}, \code{mu.lam}, \code{sig.sq.lam}
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
  \item{RJ}{
logical value indicating whether the algorithm allows the number of time partitions and 
their positions to randomly change using a reversible jump MCMC.
}
}
\details{
\tabular{ll}{
\code{y} \tab a vector of \code{n} times to the event \cr
\code{delta} \tab a vector of \code{n} censoring indicators for the event time (1=event occurred, 0=censored) \cr
\code{x} \tab covariate matrix, \code{n} observations by \code{p} variables\cr
\code{c.lam} \tab the prior parameter that controls the level of dependence between 
	adjacent intervals \cr
	\tab under Gaussian intrinsic conditional autoregression (ICAR) prior specification for \eqn{\lambda}. \eqn{0\le c.lam \le 1}\cr
\code{C} \tab the upper bound for the sum of the "birth" and "death" move probabilities for \eqn{s}, \eqn{0< C < 1}\cr
\code{a} \tab the shape parameter of conjugate gamma prior for \eqn{1/\sigma_{\lambda}^2}\cr
\code{b} \tab the rate parameter of conjugate gamma prior for \eqn{1/\sigma_{\lambda}^2}\cr
\code{alpha} \tab the rate parameter of Poisson prior for the number of time splits \eqn{J}\cr
\code{delPert} \tab the parameter to control the perturbation in the "birth" move. \eqn{0< delPert \le 0.5}\cr
\code{smax} \tab the largest observed event time \cr
\code{s.propBI} \tab a set of proposal values for the "birth" move of \eqn{s}\cr
\code{beta.ini} \tab the starting value for \eqn{\beta}\cr
\code{s} \tab the starting value for the time splits \eqn{s}\cr
\code{mu.lam} \tab the starting values for \eqn{\mu_{\lambda}}\cr
\code{sig.sq.lam} \tab the starting values for \eqn{\sigma_{\lambda}^2}\cr
}
}

\value{
\code{BpeSurv} returns an object of class \code{BpeSurv} \cr
\item{move.ind}{selected moves in the MCMC sampling}
\item{beta.p}{posterior samples for \eqn{\beta}}
\item{lambda.p}{posterior samples for \eqn{\lambda}}
\item{s.p}{posterior samples for \eqn{s}}
\item{mu.lam.p}{posterior samples for \eqn{\mu_{\lambda}}}
\item{sig.sq.lam.p}{posterior samples for \eqn{\sigma_{\lambda}^2}}
\item{accept.beta}{the number of acceptance in sampling \eqn{\beta}}
\item{accept.lambda}{the number of acceptance in sampling \eqn{\lambda}}
\item{accept.bi}{the number of acceptance in 'birth' move}
\item{accept.di}{the number of acceptance in 'death' move}
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

	priorPara$c.lam		<- 1			
	priorPara$C			<- 0.70				
	priorPara$a			<- 0.5	# prior parameter for sigma^(-2)
	priorPara$b			<- 0.01
	priorPara$alpha		<- 5	# prior parameter for J-1
	priorPara$delPert	<- 0.5   # perturbation U ~ unif(0.5 - delPert, 0.5 + delPert)
	
	initial$s				<- c(seq(8, max(y[delta == 1]), 8), max(y[delta == 1]));
		
	priorPara$smax			<- max(initial$s)					
	priorPara$s.propBI 		<- unique(sort(c(seq(0.5, priorPara$smax, 0.5),  priorPara$smax)))
					
	initial$J				<- length(initial$s)-1;
	initial$beta.ini		<- rep(1, p)
	initial$lambda			<- runif(initial$J+1, -5, 0);
	initial$sig.sq.lam		<- var(initial$lambda)
			
	initial$mu.lam			<- mean(initial$lambda)
	
	
	# fitting semi-parametric regression model to right censored survival datas	
	# set 'num.reps' to a larger number such as 500000	
	
	fit.BpeSurv <- BpeSurv(survObj, priorPara, initial, num.reps = 100,
	 						thin = 10, chain = 1, save = 100, RJ = TRUE)
	
	# check the directory '/mcmcOutcome' under your working directory
		
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ semi-parametric analysis }
\keyword{ univariate survival analysis }% __ONLY ONE__ keyword per line
