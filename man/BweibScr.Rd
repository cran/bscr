\name{BweibScr}
\alias{BweibScr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to fit a Bayesian parametric (Weibull) regression model to semi-competing risks data.
}
\description{
The function to fit a Bayesian parametric (Weibull) regression model to semi-competing risks data.
}
\usage{
BweibScr(survObj, priorPara, initial, num.reps, thin = 1, chain = 1,
			save = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{survObj}{
				The list containing observed data from \code{n} subjects; 
				\code{y1}, \code{y2}, \code{delta1}, \code{delta2}, \code{x}
				}
\item{priorPara}{
				The list containing prior parameter values; 
					\code{a1}, \code{b1}, \code{a2}, \code{b2}, \code{a3}, \code{b3},
					\code{c1}, \code{d1}, \code{c2}, \code{d2}, \code{c3}, \code{d3},					
					\code{psi}, \code{omega}
					}
					
\item{initial}{
				The list containing the starting values of the parameters;
				\code{beta1.ini}, \code{beta2.ini}, \code{beta3.ini},
				\code{alpha1.ini}, \code{alpha2.ini}, \code{alpha3.ini},
				\code{kappa1.ini}, \code{kappa2.ini}, \code{kappa3.ini},
				\code{gamma.ini}, \code{theta.ini}
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
\code{y1} \tab a vector of \code{n} times to non-terminal event \cr
\code{y2} \tab a vector of \code{n} times to terminal event \cr
\code{delta1} \tab a vector of \code{n} censoring indicators for the non-terminal event time (1=event occurred, 0=censored) \cr
\code{delta2} \tab a vector of \code{n} censoring indicators for the terminal event time (1=event occurred, 0=censored)\cr
\code{x} \tab covariate matrix, \code{n} observations by \code{p} variables\cr
\code{a1} \tab the shape parameter of conjugate gamma prior for \eqn{\alpha_1}\cr
\code{b1} \tab the rate parameter of conjugate gamma prior for \eqn{\alpha_1}\cr
\code{a2} \tab the shape parameter of conjugate gamma prior for \eqn{\alpha_2}\cr
\code{b2} \tab the rate parameter of conjugate gamma prior for \eqn{\alpha_2}\cr
\code{a3} \tab the shape parameter of conjugate gamma prior for \eqn{\alpha_3}\cr
\code{b3} \tab the rate parameter of conjugate gamma prior for \eqn{\alpha_3}\cr
\code{c1} \tab the shape parameter of conjugate gamma prior for \eqn{\kappa_1}\cr
\code{d1} \tab the rate parameter of conjugate gamma prior for \eqn{\kappa_1}\cr
\code{c2} \tab the shape parameter of conjugate gamma prior for \eqn{\kappa_2}\cr
\code{d2} \tab the rate parameter of conjugate gamma prior for \eqn{\kappa_2}\cr
\code{c3} \tab the shape parameter of conjugate gamma prior for \eqn{\kappa_3}\cr
\code{d3} \tab the rate parameter of conjugate gamma prior for \eqn{\kappa_3}\cr
\code{psi} \tab the shape parameter of the gamma prior for \eqn{1/\theta}\cr
\code{omega} \tab the rate parameter of the gamma prior for \eqn{1/\theta}\cr
\code{indGamma} \tab indices of subjects for whom the posterior samples of \eqn{\gamma} are saved\cr
\code{beta1.ini} \tab the starting values for \eqn{\beta_1}\cr
\code{beta2.ini} \tab the starting values for \eqn{\beta_2}\cr
\code{beta3.ini} \tab the starting values for \eqn{\beta_3}\cr
\code{alpha1.ini} \tab the starting values for \eqn{\alpha_1}\cr
\code{alpha2.ini} \tab the starting values for \eqn{\alpha_2}\cr
\code{alpha3.ini} \tab the starting values for \eqn{\alpha_3}\cr
\code{kappa1.ini} \tab the starting values for \eqn{\kappa_1}\cr
\code{kappa2.ini} \tab the starting values for \eqn{\kappa_2}\cr
\code{kappa3.ini} \tab the starting values for \eqn{\kappa_3}\cr
\code{gamma.ini} \tab the starting values for \eqn{\gamma}\cr
\code{theta.ini} \tab the starting values for \eqn{\theta}\cr
}
}

\value{
\code{BweibScr} returns an object of class \code{BweibScr} \cr
\item{move.ind}{selected moves in the MCMC sampling}
\item{beta1.p}{posterior samples for \eqn{\beta_1}}
\item{beta2.p}{posterior samples for \eqn{\beta_2}}
\item{beta3.p}{posterior samples for \eqn{\beta_3}}
\item{alpha1.p}{posterior samples for \eqn{\alpha_1}}
\item{alpha2.p}{posterior samples for \eqn{\alpha_2}}
\item{alpha3.p}{posterior samples for \eqn{\alpha_3}}
\item{kappa1.p}{posterior samples for \eqn{\kappa_1}}
\item{kappa2.p}{posterior samples for \eqn{\kappa_2}}
\item{kappa3.p}{posterior samples for \eqn{\kappa_3}}
\item{gamma.p}{posterior samples for \eqn{\gamma}}
\item{theta.p}{posterior samples for \eqn{\theta}}
\item{accept.beta1}{the number of acceptance in sampling \eqn{\beta_1}}
\item{accept.beta2}{the number of acceptance in sampling \eqn{\beta_2}}
\item{accept.beta3}{the number of acceptance in sampling \eqn{\beta_3}}
\item{accept.alpha1}{the number of acceptance in sampling \eqn{\alpha_1}}
\item{accept.alpha2}{the number of acceptance in sampling \eqn{\alpha_2}}
\item{accept.alpha3}{the number of acceptance in sampling \eqn{\alpha_3}}
\item{accept.theta}{the number of acceptance in sampling \eqn{\theta}}
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2013). 
Bayesian Semi-parametric Analysis of Semi-competing Risks Data: 
Estimating Readmission Rates among Pancreatic Cancer Patients, available from the authors upon request.
}
\author{
Kyu Ha Lee, Sebastien Haneuse, Deborah Schrag, and Francesca Dominici \cr
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
   				   

	# Set prior parameters and initial values
	
	survObj	<- priorPara	<- initial	<- list()
		
	y1 <- survObj$y1 <- simData$y1
	y2 <- survObj$y2 <- simData$y2	
		
	delta1 <- survObj$delta1 <- simData$delta1
	delta2 <- survObj$delta2 <- simData$delta2	
	survObj$x	<- x

	priorPara$a1			<- 0.5	# prior parameters for alpha
	priorPara$b1			<- 0.01
	priorPara$a2			<- 0.5
	priorPara$b2			<- 0.01
	priorPara$a3			<- 0.5
	priorPara$b3			<- 0.01
	
	priorPara$c1			<- 0.5	# prior parameters for kappa
	priorPara$d1			<- 0.05
	priorPara$c2			<- 0.5
	priorPara$d2			<- 0.05
	priorPara$c3			<- 0.5
	priorPara$d3			<- 0.05
				
	priorPara$psi		<- 0.7		# prior parameter for theta
	priorPara$omega		<- 0.7		
	
	priorPara$indGamma <- 1:10		
		
	initial$beta1.ini		<- rep(1, p)
	initial$beta2.ini		<- rep(1, p)
	initial$beta3.ini		<- rep(1, p)
		
	initial$alpha1.ini			<- runif(1, 0, 2);
	initial$alpha2.ini			<- runif(1, 0, 2);
	initial$alpha3.ini			<- runif(1, 0, 2);
		
	initial$kappa1.ini			<- 1
	initial$kappa2.ini			<- 1
	initial$kappa3.ini			<- 1
					
	initial$theta.ini		<- 1
	
	initial$gamma.ini		<- rgamma(n, shape = 1/initial$theta.ini, rate = 1/initial$theta.ini);
	
		
	# fitting parametric regression model to semi-competing risks data	
	# set 'num.reps' to a larger number such as 500000	
	
	fit.BweibScr <- BweibScr(survObj, priorPara, initial, num.reps = 100,
							thin = 10, chain = 1, save = 100)
	
		
	# check the directory '/mcmcOutcome' under your working directory
		
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ parametric analysis }
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line
