\name{EhrBweib}
\alias{EhrBweib}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to calculate the conditional explanatory hazard ratio (EHR) under parametric model
}
\description{
The function to calculate the conditional explanatory hazard ratio (EHR) under parametric model
}
\usage{
EhrBweib(x, beta2.p, beta3.p, alpha2.p, alpha3.p, kappa2.p, kappa3.p, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a covariate vector of length \code{p} for a subject
}
  \item{beta2.p}{
posterior samples for \eqn{\beta_2}
}
  \item{beta3.p}{
posterior samples for \eqn{\beta_3}
}
  \item{alpha2.p}{
posterior samples for \eqn{\alpha_2}
}
  \item{alpha3.p}{
posterior samples for \eqn{\alpha_3}
}
  \item{kappa2.p}{
posterior samples for \eqn{\kappa_2}
}
  \item{kappa3.p}{
posterior samples for \eqn{\kappa_3}
}
  \item{s}{
time points where EHR is computed
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a matrix of calculated values where each row corresponds to a value of conditional EHR
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2013). 
Bayesian Semi-parametric Analysis of Semi-competing Risks Data: 
Estimating Readmission Rates among Pancreatic Cancer Patients, available from the authors upon request.
}
\author{
Kyu Ha Lee, Sebastein Haneuse, Deborah Schrag, and Francesca Dominici \cr
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
	
	fit.BweibScr <- BweibScr(survObj, priorPara, initial, num.reps = 500,
							thin = 10, chain = 1, save = 500)
							
	
	beta2.p		<- fit.BweibScr$beta2.p
	beta3.p		<- fit.BweibScr$beta3.p
	alpha2.p	<- fit.BweibScr$alpha2.p
	alpha3.p	<- fit.BweibScr$alpha3.p
	kappa2.p	<- fit.BweibScr$kappa2.p
	kappa3.p	<- fit.BweibScr$kappa3.p	

	s 		<- 1:40


	ehr <- EhrBweib(x[2,], beta2.p, beta3.p, alpha2.p, alpha3.p, kappa2.p, kappa3.p, s)
	
	ehrMed		<- apply(ehr, 2, median)
	ehr0.975	<- apply(ehr, 2, quantile, prob = 0.975)
	ehr0.025	<- apply(ehr, 2, quantile, prob = 0.025)


	plot(s, ehrMed, type = "l", ylim = c(0, 5))
	lines(s, ehr0.975, lty = 2)
	lines(s, ehr0.025, lty = 2)		
							
							
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ parametric analysis }
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line

