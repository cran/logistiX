\name{print.logistiX}
\alias{print.logistiX}
\title{
Print method for \code{logistiX}
}
\description{
This method prints the main results from a \code{logistiX} analysis.

}
\usage{
\method{print}{logistiX}(x,...)
}
\arguments{
  \item{x}{
  a \code{logistiX} object.
}
\item{...}{further arguments}
}
\value{
returns a 4p x 5 matrix of estimates and confidence limits, where p is the number of variables. For each variable, four methods are applied:

for estimation: MLE, MUE, LX, and CCFL; for confidence limits: TST, TST-pmid, scores, scores-pmid
}
\references{
Heinze G and Puhr R (2010). Bias-reduced and separation-proof conditional logistic regression with small or sparse data sets. Statistics in Medicine 29:770-777
}
\author{
Georg Heinze and Tobias Ladner
}
\note{
\code{\link{summary.logistiX}} provides a nicer summary with user-specifiable contents.
}


\seealso{
\code{\link{summary.logistiX}}
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

  print(test1)
}  
\keyword{regression}
\keyword{models}
