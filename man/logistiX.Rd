\name{logistiX}
\alias{logistiX}
\title{
Exact logistic regression
}
\description{
Main function for exact logistic regression in R.
}
\usage{
logistiX(x = NULL, y = NULL, strat = NULL, tst = 1, sc = 1, pmid = 2, option = 2, noint = 0, vstart = 1, vstop = NA, alpha = 0.05, details = 0)
}

\arguments{
  \item{x}{
A n x p matrix of covariates, these must be 0/1 coded. (Currently, it is not possible to use discrete/ordinal variables or continuous variables.)
}
  \item{y}{
A p x 1 vector of binary outcomes
}
  \item{strat}{
A p x 1 vector of stratification categories (e.g., matched set indicator in a matched case-control study). Not needed.
}
  \item{tst}{
Compute twice-smaller-tail method of inference (default=TRUE=yes)
}
  \item{sc}{
Compute scores method of inference (default=TRUE=yes)
}
  \item{pmid}{
Request pmid adjustment (default=2 meaning that pmid and exact results are computed. 1: only pmid, 0: only exact)
}
  \item{option}{
 option=2: compute permutational distributions and estimate; option=1 only compute permutational distributions
}
  \item{noint}{
 Omit intercept (NOT WORKING!) (default=0)
 }
  \item{vstart}{
First variable to estimate, default=1
}
  \item{vstop}{
Last variable to estimate, default = p
}
  \item{alpha}{
1-confidence level, default=0.05
}
  \item{details}{
  Request detailed printing of messages (default=0=no).
}
}
\details{

\code{logistiX} computes the permutational null distribution of the sufficient statistics of either all or a subset of the variables of a model, using the multivariate shift algorithm (Hirji et al 1987).
By default, \code{logistiX} also computes the parameter estimates by four methods: MLE (maximum conditional likelihood), MUE (median unbiased estimation), 
LX (the "LogXact" method, which uses MLE if defined, MUE otherwise), CCFL (Firth-corrected estimate, Heinze and Puhr, 2010). Confidence intervals are computed by the following methods:
TST: twice-smaller-tail method, scores: scores method, pmid: adds pmid-adjustment (the observed sufficient statistic is downweighted by half (the "pmid-factor"), a correction to 
remove conservatism by exact inference; these confidence intervals are 
only approximately exact), randomized: computes integrated randomized confidence intervals, i.e., confidence limits are averaged over a series of values obtained
by varying the pmid-factor from 0 to 1.

By calling the summary, coef or confint methods, one can obtain results by the specified methods.
}
\value{
  A list with components:
  \item{distout}{Conditional permutational distribution for variables \code{vstart:vstop}}
  \item{estout}{Parameter estimates by the methods MLE, MUE, LX and CCFL}
  \item{ciout}{Confidence intervals by the methods TST, TST-pmid, SC and SC-pmid}
  \item{tobs}{Observed values of the sufficient statistics for all variables (including intercept)}
  \item{call}{The call to \code{logistiX}}
%% ...
}
\references{
Hirji KF, Mehta CR and Patel NR (1987). Computing Distributions for Exact Logistic Regression. Journal of the American Statistical Association 82:1110-1117

Heinze G and Puhr R (2010). Bias-reduced and separation-proof conditional logistic regression with small or sparse data sets. Statistics in Medicine 29:770-777
}
\author{
Georg Heinze and Tobias Ladner 
(georg.heinze@meduniwien.ac.at)
}
\note{
Currently, the implementation has two major limitations:
  (1) only binary covariates can be processed, and (2) it is not possible to compute the joint distribution of several parameters and consequently no simultaneous testing is possible 
  }

\section{Warning}{ 
  Memory requirements can be excessive already with a moderately-sized data set (e.g., with N=200 and 5 variables).
  R's 64-bit version allows for more allocatable memory, but can be slow if virtual memory is on a conventional harddisk drive.
  In such cases, the \code{elrm} package may be an alternative.
  }

\seealso{
\code{\link{print.logistiX}},
\code{\link{summary.logistiX}},
\code{\link{plot.logistiX}},
\code{\link{profile.logistiX}},
\code{\link{coef.logistiX}},
\code{\link{confint.logistiX}}
 }
\examples{
 set.seed(123)
  vars<-5
  n<-50
  beta<-c(-2,3,2,-2,1,2,2,1,0,-1,2)
  px<-0.2

  x<-matrix(rbinom(vars*n,1,px),n,vars)

  py<-1/(1+exp(-cbind(rep(1,n),x) \%*\% beta[1:(vars+1)]))
  y<-rbinom(n,1,py)
  
  test1<-logistiX(x=x, y=y)
  summary(test1)
  
  #LogXact estimate:
  coef(test1, "LX")
  confint(test1, type="exact")

  
  # CCFL estimate:
  coef(test1, "CCFL")
  # PCCL confidence limits:
  confint(test1, type="PCCL")
}  
\keyword{regression}
\keyword{models}
