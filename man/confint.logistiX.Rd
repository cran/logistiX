\name{confint.logistiX}
\alias{confint.logistiX}
\title{
Method to extract confidence intervals from a \code{logistiX} object.   }
\description{
This method extracts confidence intervals from a \code{logistiX} object. 

}
\usage{
\method{confint}{logistiX}(object, parm, level = 0.95,  type = "exact", method = "TST", ran.steps = seq(5/100, 0.95, 0.01), interval = c(-19, 19), vstart = NULL, vstop = NULL, ...)
}
\arguments{
  \item{object}{
 a \code{logistiX} object}
  \item{parm}{
  a specification of which parameters are to be given confidence intervals,  a vector of numbers. If missing, all parameters are considered. \code{vstart} and \code{vstop} may be used alternatively.)
  }
   \item{level}{
  the confidence level (default=0.95)
  }

  \item{type}{
   one of \code{c("exact", "pmid","randomized","PCCL")} (see above)%%     ~~Describe \code{type} here~~
}
  \item{method}{
  one of \code{c("TST", "scores")} (see above)
}
  \item{ran.steps}{
  the pmid-factor steps used for integration (\code{randomized} CL type). Default is 0.05 to 0.95 by 0.01 (\code{ran.steps = seq(5/100, 0.95, 0.01)}).
}
  \item{interval}{
  the range of values for \code{optimize} to search the confidence limits in. If estimation fails (i.e., if the actual level is not equal the specified level), then the interval limits
  are halfened.
  

}
  \item{vstart}{
  Optional: the first variable to compute confidence intervals (default=\code{NULL}). If \code{NULL}, the first variable contained in the \code{obj} object will be used.
  }
  \item{vstop}{
  Optional: the last variable to compute confidence intervals (default=\code{NULL}).  If \code{NULL}, the last variable contained in the \code{obj} object will be used.
}
\item{...}{
additional argument(s) for methods
}
}
\details{
This \code{confint} method can be used to generate confidence intervals by various  methods and 'types', that can be combined. Available methods are:

\tabular{lll}{
Method \tab \code{type} argument \tab \code{method} argument\cr
Exact twice-smaller-tail \tab \code{"exact"} \tab \code{"TST"} \cr
Pmid twice-smaller-tail \tab \code{"pmid"} \tab \code{"TST"}  \cr
Exact scores \tab \code{"exact"} \tab \code{"scores"} \cr
Pmid scores \tab \code{"pmid"} \tab \code{"scores"}  \cr
Integrated randomized twice-smaller-tail (EXPERIMENTAL)\tab \code{"randomized"} \tab \code{"TST"} \cr
Profile penalized completely conditional likelihood \tab \code{"PCCL"} \tab \cr
}

The confidence level can be chosen, and for the \code{"randomized"} method, the number and density of integration steps can be specified. Furthermore,
the set of parameters for which confidence limits are requested can be restricted.

PCCL intervals are as proposed by Heinze and Puhr (2010).
Integrated randomized TST intervals are an experimental feature. 

The \code{confint} method is also called by \code{summary}. 
}
\value{
A p x 2 matrix of confidence limits, where p = \code{vstop}-\code{vstart}+1.
}
\references{
Heinze G and Puhr R (2010). Bias-reduced and separation-proof conditional logistic regression with small or sparse data sets. Statistics in Medicine 29:770-777
}
\author{
Georg Heinze and Tobias Ladner
}
\note{
Comparison to LogXact: in LogXact, only exact TST and PCCL intervals are implemented.
}


\seealso{
\code{\link{coef.logistiX}}
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
  confint(test1, type="exact")
  confint(test1, type="pmid")
}
\keyword{regression}
\keyword{models}
