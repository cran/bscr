\name{bscr}
\alias{bscr}
\docType{package}
\title{
Algorithms for fitting Bayesian parametric and semi-parametric models to semi-competing risks data.
}
\description{
The package provides Bayesian schemes for fitting the following survival models: \cr \cr
1) Bayesian semi-parametric model for semi-competing risks data, \cr
2) Bayesian parametric model for semi-competing risks data, \cr
3) Bayesian semi-parametric model for univariate right censored survival data,\cr
4) Bayesian parametric model for univariate right censored survival data \cr
}
\details{
The package includes following functions:
\tabular{ll}{
\code{BpeScr} \tab The function to fit Bayesian piecewise exponential (semi-parametric) model \cr
				\tab to semi-competing risks data \cr
\code{BweibScr} \tab The function to fit Bayesian weibull (parametric) model to semi-competing risks data \cr
\code{BpeSurv} \tab  The function to fit Bayesian Bayesian piecewise exponential (semi-parametric) \cr 
				\tab to univariate right censored survival data \cr
\code{BweibSurv} \tab The function to fit Bayesian weibull (parametric) model to univariate right censored survival data \cr
\code{EhrBpe} \tab  The function to calculate the conditional explanatory hazard ratio (EHR) \cr
				\tab under semi-parametric model \cr
\code{EhrBweib} \tab  The function to calculate the conditional explanatory hazard ratio (EHR) \cr
				\tab under parametric model \cr
\code{SimScr} \tab The function to simulate semi-competing risks data under Weibull model \cr
\code{SimSurv} \tab The function to simulate right censored survival data under Weibull model \cr
}
\tabular{ll}{
Package: \tab bscr\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-01-02\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Kyu Ha Lee, Sebastein Haneuse, Deborah Schrag, and Francesca Dominici \cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
\references{
Lee, K. H., Haneuse, S., Schrag, D., and Dominici, F. (2013). 
Bayesian Semi-parametric Analysis of Semi-competing Risks Data: 
Estimating Readmission Rates among Pancreatic Cancer Patients, available from the authors upon request.
}
%Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%documentation directory ~~
\keyword{ package }
%\seealso{
%Optional links to other man pages, e.g. ~~
%\code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%
%}
