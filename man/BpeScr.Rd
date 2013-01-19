\name{BpeScr}
\alias{BpeScr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to fit a Bayesian semi-parametric regression model to semi-competing risks data.
}
\description{
The function to fit a Bayesian semi-parametric regression model to semi-competing risks data.
}
\usage{
BpeScr(survObj, priorPara, initial, num.reps, thin = 1, chain = 1,
		save = 1000, RJ = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{survObj}{
				The list containing observed data from \code{n} subjects; 
				\code{y1}, \code{y2}, \code{delta1}, \code{delta2}, \code{x}
				}
\item{priorPara}{
				The list containing prior parameter values; 
					\code{c.lam1}, \code{c.lam2}, \code{c.lam3}, \code{C1}, \code{C2}, \code{C3},
					\code{a1}, \code{b1}, \code{a2}, \code{b2}, \code{a3}, \code{b3},
					\code{alpha1}, \code{alpha2}, \code{alpha3}, \code{psi}, \code{omega},
					\code{delPert1}, \code{delPert2}, \code{delPert3}, \code{s1max},
					\code{s2max}, \code{s3max}, \code{s.propBI1}, \code{s.propBI2},
					, \code{s.propBI3}
					}
					
\item{initial}{
				The list containing the starting values of the parameters;
				\code{beta1.ini}, \code{beta2.ini}, \code{beta3.ini},
				\code{s1}, \code{s2}, \code{s3}, \code{lambda1.ini}, \code{lambda2.ini}, \code{lambda3.ini},
				\code{mu.lam1}, \code{sig.sq.lam1}, \code{mu.lam2}, \code{sig.sq.lam2}, \code{mu.lam3}, \code{sig.sq.lam3},
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
  \item{RJ}{
logical value indicating whether the algorithm allows the number of time partitions and 
their positions to randomly change using a reversible jump MCMC.
}
}
\details{
\tabular{ll}{
\code{y1} \tab a vector of \code{n} times to non-terminal event \cr
\code{y2} \tab a vector of \code{n} times to terminal event \cr
\code{delta1} \tab a vector of \code{n} censoring indicators for the non-terminal event time (1=event occurred, 0=censored) \cr
\code{delta2} \tab a vector of \code{n} censoring indicators for the terminal event time (1=event occurred, 0=censored)\cr
\code{x} \tab covariate matrix, \code{n} observations by \code{p} variables\cr
\code{c.lam1} \tab the prior parameter that controls the level of dependence between 
	adjacent intervals \cr
	\tab under Gaussian intrinsic conditional autoregression (ICAR) prior specification for \eqn{\lambda_1}.\cr
	\tab \eqn{0\le c.lam1 \le 1}\cr
\code{c.lam2} \tab the prior parameter that controls the level of dependence between 
	adjacent intervals \cr
	\tab under Gaussian intrinsic conditional autoregression (ICAR) prior specification for \eqn{\lambda_2}.\cr
	\tab \eqn{0\le c.lam2 \le 1}\cr
\code{c.lam3} \tab the prior parameter that controls the level of dependence between 
	adjacent intervals \cr
	\tab under Gaussian intrinsic conditional autoregression (ICAR) prior specification for \eqn{\lambda_3}.\cr
	\tab \eqn{0\le c.lam3 \le 1}\cr
			
\code{C1} \tab the upper bound for the sum of the "birth" and "death" move probabilities for \eqn{s_1}\cr
\code{C2} \tab the upper bound for the sum of the "birth" and "death" move probabilities for \eqn{s_2}\cr
\code{C3} \tab the upper bound for the sum of the "birth" and "death" move probabilities for \eqn{s_3}\cr
			\tab \eqn{0 < C1+C2+C3 < 1}\cr
\code{a1} \tab the shape parameters of conjugate gamma prior for \eqn{1/\sigma_{\lambda_1}^2}\cr
\code{b1} \tab the rate parameters of conjugate gamma prior for \eqn{1/\sigma_{\lambda_1}^2}\cr
\code{a2} \tab the shape parameters of conjugate gamma prior for \eqn{1/\sigma_{\lambda_2}^2}\cr
\code{b2} \tab the rate parameters of conjugate gamma prior for \eqn{1/\sigma_{\lambda_2}^2}\cr
\code{a3} \tab the shape parameters of conjugate gamma prior for \eqn{1/\sigma_{\lambda_3}^2}\cr
\code{b3} \tab the rate parameters of conjugate gamma prior for \eqn{1/\sigma_{\lambda_3}^2}\cr
\code{alpha1} \tab the rate parameter of Poisson prior for the number of time splits \eqn{J_1}\cr
\code{alpha2} \tab the rate parameter of Poisson prior for the number of time splits \eqn{J_2}\cr
\code{alpha3} \tab the rate parameter of Poisson prior for the number of time splits \eqn{J_3}\cr
\code{psi} \tab the shape parameter of the gamma prior for \eqn{1/\theta}\cr
\code{omega} \tab the rate parameter of the gamma prior for \eqn{1/\theta}\cr
\code{delPert1} \tab the parameter to control the perturbation in the "birth" move for \eqn{s_1}. \eqn{0< delPert1 \le 0.5}\cr
\code{delPert2} \tab the parameter to control the perturbation in the "birth" move for \eqn{s_2}. \eqn{0< delPert1 \le 0.5}\cr
\code{delPert3} \tab the parameter to control the perturbation in the "birth" move for \eqn{s_3}. \eqn{0< delPert1 \le 0.5}\cr
\code{indGamma} \tab indices of subjects for whom the posterior samples of \eqn{\gamma} are saved\cr
\code{s1max} \tab the largest observed event times for transition \eqn{g = 1}\cr
\code{s2max} \tab the largest observed event times for transition \eqn{g = 2}\cr
\code{s3max} \tab the largest observed event times for transition \eqn{g = 3}\cr
\code{s.propBI1} \tab a set of proposal values for the "birth" move of \eqn{s_1}\cr
\code{s.propBI2} \tab a set of proposal values for the "birth" move of \eqn{s_2}\cr
\code{s.propBI3} \tab a set of proposal values for the "birth" move of \eqn{s_3}\cr
\code{beta1.ini} \tab the starting values for \eqn{\beta_1}\cr
\code{beta2.ini} \tab the starting values for \eqn{\beta_2}\cr
\code{beta3.ini} \tab the starting values for \eqn{\beta_3}\cr
\code{s1} \tab the starting values for the time splits \eqn{s_1}\cr
\code{s2} \tab the starting values for the time splits \eqn{s_2}\cr
\code{s3} \tab the starting values for the time splits \eqn{s_3}\cr
\code{mu.lam1} \tab the starting values for \eqn{\mu_{\lambda_1}}\cr
\code{sig.sq.lam1} \tab the starting values for \eqn{\sigma_{\lambda_1}^2}\cr
\code{mu.lam2} \tab the starting values for \eqn{\mu_{\lambda_2}}\cr
\code{sig.sq.lam2} \tab the starting values for \eqn{\sigma_{\lambda_2}^2}\cr
\code{mu.lam3} \tab the starting values for \eqn{\mu_{\lambda_3}}\cr
\code{sig.sq.lam3} \tab the starting values for \eqn{\sigma_{\lambda_3}^2}\cr
\code{gamma.ini} \tab the starting values for \eqn{\gamma}\cr
\code{theta.ini} \tab the starting values for \eqn{\theta}\cr

}
}

\value{
\code{BpeScr} returns an object of class \code{BpeScr} \cr
\item{move.ind}{selected moves in the MCMC sampling}
\item{beta1.p}{posterior samples for \eqn{\beta_1}}
\item{beta2.p}{posterior samples for \eqn{\beta_2}}
\item{beta3.p}{posterior samples for \eqn{\beta_3}}
\item{lambda1.p}{posterior samples for \eqn{\lambda_1}}
\item{lambda2.p}{posterior samples for \eqn{\lambda_2}}
\item{lambda3.p}{posterior samples for \eqn{\lambda_3}}
\item{s1.p}{posterior samples for \eqn{s_1}}
\item{s2.p}{posterior samples for \eqn{s_2}}
\item{s3.p}{posterior samples for \eqn{s_3}}
\item{mu.lam1.p}{posterior samples for \eqn{\mu_{\lambda_1}}}
\item{mu.lam2.p}{posterior samples for \eqn{\mu_{\lambda_2}}}
\item{mu.lam3.p}{posterior samples for \eqn{\mu_{\lambda_3}}}
\item{sig.sq.lam1.p}{posterior samples for \eqn{\sigma_{\lambda_1}^2}}
\item{sig.sq.lam2.p}{posterior samples for \eqn{\sigma_{\lambda_2}^2}}
\item{sig.sq.lam3.p}{posterior samples for \eqn{\sigma_{\lambda_3}^2}}
\item{gamma.p}{posterior samples for \eqn{\gamma}}
\item{theta.p}{posterior samples for \eqn{\theta}}
\item{accept.beta1}{the number of acceptance in sampling \eqn{\beta_1}}
\item{accept.beta2}{the number of acceptance in sampling \eqn{\beta_2}}
\item{accept.beta3}{the number of acceptance in sampling \eqn{\beta_3}}
\item{accept.lambda1}{the number of acceptance in sampling \eqn{\lambda_1}}
\item{accept.lambda2}{the number of acceptance in sampling \eqn{\lambda_2}}
\item{accept.lambda3}{the number of acceptance in sampling \eqn{\lambda_3}}
\item{accept.theta}{the number of acceptance in sampling \eqn{\theta}}
\item{accept.bi1}{the number of acceptance in 'birth' moves for for \eqn{s_1}}
\item{accept.bi2}{the number of acceptance in 'birth' moves for for \eqn{s_2}}
\item{accept.bi3}{the number of acceptance in 'birth' moves for for \eqn{s_3}}
\item{accept.di1}{the number of acceptance in 'death' moves for for \eqn{s_1}}
\item{accept.di2}{the number of acceptance in 'death' moves for for \eqn{s_2}}
\item{accept.di3}{the number of acceptance in 'death' moves for for \eqn{s_3}}
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
	
	fit.BpeScr <- BpeScr(survObj, priorPara, initial, num.reps = 100, thin = 10, 
						chain = 1, save = 100, RJ = TRUE)
	
	# check the directory '/mcmcOutcome' under your working directory
	

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ semi-parametric analysis }
\keyword{ semi-competing risks analysis }% __ONLY ONE__ keyword per line
