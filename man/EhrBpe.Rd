\name{EhrBpe}
\alias{EhrBpe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to calculate the conditional explanatory hazard ratio (EHR) under semi-parametric model
}
\description{
The function to calculate the conditional explanatory hazard ratio (EHR) under semi-parametric model
}
\usage{
EhrBpe(x, beta2.p, beta3.p, lambda2.p, lambda3.p, s2.p, s3.p, s)
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
  \item{lambda2.p}{
posterior samples for \eqn{\lambda_2}
}
  \item{lambda3.p}{
posterior samples for \eqn{\lambda_3}
}
  \item{s2.p}{
posterior samples for \eqn{s_2}
}
  \item{s3.p}{
posterior samples for \eqn{s_3}
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
	alpha3.true = 0.5
		
	kappa1.true = 0.02
	kappa2.true = 0.1
	kappa3.true = 0.15
	
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

	priorPara$c.lam1		<- 1
	priorPara$c.lam2		<- 1
	priorPara$c.lam3		<- 1
			
	priorPara$C1			<- 0.20
	priorPara$C2			<- 0.20
	priorPara$C3			<- 0.20
			
	priorPara$a1			<- 0.5	# prior parameter for sigma^(-2)
	priorPara$b1			<- 0.01
	priorPara$a2			<- 0.5
	priorPara$b2			<- 0.01
	priorPara$a3			<- 0.5
	priorPara$b3			<- 0.01
	
	priorPara$alpha1		<- 5	# prior parameter for J-1
	priorPara$alpha2		<- 5
	priorPara$alpha3		<- 5
			
	priorPara$psi		<- 0.7		# prior parameter for theta
	priorPara$omega		<- 0.7
	
	priorPara$delPert1	<- 0.5   # perturbation U ~ unif(0.5 - delPert, 0.5 + delPert)
	priorPara$delPert2	<- 0.5
	priorPara$delPert3	<- 0.5		
	
	priorPara$indGamma <- 1:10	
	
	initial$s1	<- c(seq(8, max(y1[delta1 == 1]), 8), max(y1[delta1 == 1]));
	initial$s2	<- c(seq(8, max(y2[delta1 == 0 & max(y1[delta2 == 1])]), 8), 
						max(y2[delta1 == 0 & max(y1[delta2 == 1])]));
	initial$s3	<- c(seq(8, max(y2[delta1 == 1 & max(y1[delta2 == 1])]), 8), 
						max(y2[delta1 == 1 & max(y1[delta2 == 1])]));
		
	priorPara$s1max			<- max(initial$s1)
	priorPara$s2max			<- max(initial$s2)
	priorPara$s3max			<- max(initial$s3)	
					
	priorPara$s.propBI1 <- unique(sort(c(seq(0.5, priorPara$s1max, 0.5),  
									priorPara$s1max)))
	priorPara$s.propBI2 <- unique(sort(c(seq(0.5, priorPara$s2max, 0.5),  
									priorPara$s2max)))
	priorPara$s.propBI3 <- unique(sort(c(seq(0.5, priorPara$s3max, 0.5), 
									priorPara$s3max)))
					
	initial$J1				<- length(initial$s1)-1;
	initial$J2				<- length(initial$s2)-1;
	initial$J3				<- length(initial$s3)-1;
		
	initial$beta1.ini		<- rep(0.5, p)
	initial$beta2.ini		<- rep(0.5, p)
	initial$beta3.ini		<- rep(0.5, p)
		
	initial$lambda1			<- runif(initial$J1+1, -5, 0);
	initial$lambda2			<- runif(initial$J2+1, -5, 0);
	initial$lambda3			<- runif(initial$J3+1, -5, 0);
						
	initial$theta.ini		<- 1	
			
	initial$sig.sq.lam1		<- var(initial$lambda1)
	initial$sig.sq.lam2		<- var(initial$lambda2)
	initial$sig.sq.lam3		<- var(initial$lambda3)
			
	initial$mu.lam1			<- mean(initial$lambda1)
	initial$mu.lam2			<- mean(initial$lambda2)
	initial$mu.lam3			<- mean(initial$lambda3)
	
	initial$gamma.ini		<- rgamma(n, shape = 1/initial$theta.ini, 
											rate = 1/initial$theta.ini);
	


	# fitting semi-parametric regression model to semi-competing risks data	
	# set 'num.reps' to a larger number such as 500000
	
	fit.BpeScr <- BpeScr(survObj, priorPara, initial, num.reps = 500, thin = 10, 
						chain = 1, save = 500, RJ = TRUE)

	beta2.p		<- fit.BpeScr$beta2.p
	beta3.p		<- fit.BpeScr$beta3.p
	lambda2.p	<- fit.BpeScr$lambda2.p
	lambda3.p	<- fit.BpeScr$lambda3.p
	s2.p		<- fit.BpeScr$s2.p
	s3.p		<- fit.BpeScr$s3.p

	s 		<- 1:min(fit.BpeScr$priorPara$s2max, fit.BpeScr$priorPara$s3max)


	ehr <- EhrBpe(x[2,], beta2.p, beta3.p, lambda2.p, lambda3.p, s2.p, s3.p, s)
	
	ehrMed		<- apply(ehr, 2, median)
	ehr0.975	<- apply(ehr, 2, quantile, prob = 0.975)
	ehr0.025	<- apply(ehr, 2, quantile, prob = 0.025)


	plot(s, ehrMed, type = "l", ylim = c(0, 7))
	lines(s, ehr0.975, lty = 2)
	lines(s, ehr0.025, lty = 2)	

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ semi-parametric analysis }
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line
