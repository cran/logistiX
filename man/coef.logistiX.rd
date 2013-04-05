\name{coef.logistiX}
\alias{coef.logistiX}
\title{
Method to extract regression coefficients from a \code{logistiX} object
}
\description{
 This method returns the estimated regression coefficients from a \code{logistiX} object. Available estimation methods are maximum conditional likelihood,
 median unbiased, their combination (median unbiased replaces infinite maximum likelihood estimates), and completely conditional Firth-type likelihood maximization (CCFL, Heinze and Puhr, 2010)}
\usage{
\method{coef}{logistiX}(object, type = "LX", ...)
}
\arguments{
  \item{object}{ 
  a \code{logistiX} object
%%     ~~Describe \code{obj} here~~
}
  \item{type}{
  one of \code{c("LX", "CCFL", "MLE", "MUE")}. \code{"MLE"} returns the maximum likelihood estimates, \code{"MUE"} the median unbiased estimates,
   \code{"LX"} returns the MLE, and replaces infinite values by the MUE, and \code{"CCFL"} returns the Firth-corrected (CCFL) estimate.
}
\item{...}{
  additional argument(s) for methods.
}
}
\details{
The parameter for the intercept term will not be computed and can not be extracted.
}
\value{
 A px1 vector of regression coefficients. Currently, the intercept estimate is not computed.
}
\references{
Heinze G and Puhr R (2010). Bias-reduced and separation-proof conditional logistic regression with small or sparse data sets. Statistics in Medicine 29:770-777
}
\author{
Georg Heinze and Tobias Ladner
}
\note{
Comparison to LogXact: basically, the same functionality is provided. In \code{logistiX}, the median unbiased estimate can also be displayed if the maximum conditional
likelihood estimate is finite.
}
\seealso{
\code{confint.logistiX}
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
  
  #LogXact estimate:
  coef(test1, "LX")
  #CCFL estimate:
  coef(test1, "CCFL")
 }
\keyword{regression}
\keyword{models}
