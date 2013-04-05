\name{summary.logistiX}
\alias{summary.logistiX}
\title{
Summary method for \code{logistiX}
}
\description{
This methods provides a table of parameter estimates, confidence intervals, and p-values from a \code{logistiX} analysis.

}
\usage{
\method{summary}{logistiX}(object, esttype = "LX", citype = "exact", cimethod = NULL, pmid = FALSE, testtype = "TST", cilevel=0.95, ...)
}
\arguments{
  \item{object}{
 A \code{logistiX} object.
 }
  \item{esttype}{
Type of estimate: see documentation of \code{\link{coef.logistiX}}}
  \item{citype}{
Type of confidence interval: see documentation of \code{\link{confint.logistiX}}}
  \item{cimethod}{
Method for confidence interval estimation: see documentation of \code{confint.logistiX}
}
  \item{pmid}{
 pmid adjustment for testing (if \code{TRUE})
}
  \item{testtype}{
Type of test (one of \code{"TST", "scores", "probability", "PCLR"}; default=\code{"TST"}):

\code{"TST"} requests twice-smaller-tail test

\code{"scores"} requests scores test

\code{"probability"} requests sum-of-smaller-probabilities test

\code{"PCLR"} requests penalized completely conditional likelihood ratio test (not exact)
}
\item{cilevel}{the confidence level (default=0.95)
}
\item{...}{further arguments}
}

\details{
The \code{summary} method allows to provide a comprehensive table of estimates, confidence intervals and p-values, where the methods used for these statistics 
can be  chosen by the user. 
}
\value{
Returns some summary of the function call, the requested methods and a p x 6 matrix with (p is the number of variables):
\item{column 1}{ estimates (depending on \code{esttype})}
\item{columns 2-3}{ confidence limits (depending on \code{citype} and \code{cimethod})}
\item{column 4}{ test statistics (depending on \code{testtype})}
\item{column 5}{ p-values (depending on \code{testtype})}
\item{column 6}{ cardinality of the exact distribution of sufficient statistic for the respective parameter}
}
\references{
Heinze G and Puhr R (2010). Bias-reduced and separation-proof conditional logistic regression with small or sparse data sets. Statistics in Medicine 29:770-777
}
\author{
Georg Heinze and Tobias Ladner
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

  summary(test1, esttype="CCFL", citype="PCCL")
}
\keyword{regression}
\keyword{models}
