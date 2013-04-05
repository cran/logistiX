\name{logistiX-package}
\alias{logistiX-package}
\docType{package}
\title{
Exact logistic regression
}
\description{
Implements exact conditional logistic regression in R.
}
\details{
\tabular{ll}{
Package: \tab logistiX\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-09-28\cr
License: \tab GPL\cr
}

\code{logistiX} implements exact conditional logistic regression in R, including maximum conditional likelihood, median unbiased and Firth-corrected estimation methods,
and twice-smaller-tail, probability and scores methods for hypothesis tests, with optional P-mid adjustment. For confidence interval estimation, various methods are available, 
such as exact, P-mid and integrated randomized variants, as well as the profile penalized completely conditional likelihood method discussed in Heinze and Puhr, 2010.
The construction of the exact conditional permutational distributions is based on the Multivariate Shift Algorithm by Hirji et al, 1989. In contrast to other, Monte-Carlo-
or Markov-chain-Monte-Carlo approaches, this is an exact method, and does not depend on Monte Carlo sample size or on convergence of an MCMC chain. 


Implemented Methods:

\tabular{ll}{
Algorithm to construct permutational distribution: \tab Multivariate Shift Algorithm\cr
Estimation methods: \tab Maximum conditional likelihood estimation (MCLE)\cr
 \tab Median unbiased estimation (MUE)\cr
 \tab "LogXact" (LX, combination of MCLE and MUE)\cr
 \tab Firth-corrected estimation (FC)\cr
Confidence interval methods: \tab Exact twice-smaller-tail (\code{TST})\cr
 \tab Exact scores (inverting the scores test) (\code{scores})\cr
 \tab Pmid adjustment  (\code{pmid})\cr
 \tab Integrated randomized confidence limits (\code{randomized})\cr
 \tab Profile penalized completely conditional likelihood (\code{PCCL})\cr
Inference methods: \tab Exact twice-smaller-tail (\code{TST})\cr
 \tab Exact scores (\code{scores})\cr
 \tab Exact probability (\code{probability})\cr 
 \tab P-mid adjustment (\code{pmid})\cr
}

}
\author{
Georg Heinze and Tobias Ladner

Maintainer: Georg Heinze <georg.heinze@meduniwien.ac.at>
}
\references{
Hirji KF, Mehta CR and Patel NR (1987). Computing Distributions for Exact Logistic Regression. Journal of the American Statistical Association 82:1110-1117

Heinze G and Puhr R (2010). Bias-reduced and separation-proof conditional logistic regression with small or sparse data sets. Statistics in Medicine 29:770-777
}
\keyword{regression}
\keyword{models}
\keyword{package}
\seealso{
Other packages for small-sample logistic regression are:

\code{elrm}  (Approximate exact inference in logistic regression using MCMC)

\code{logistf}  (Bias-reduced logistic regression using Firth's correction)
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
  coef(test1)
  confint(test1)

}
