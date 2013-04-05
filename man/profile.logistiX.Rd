\name{profile.logistiX}
\alias{profile.logistiX}
\title{
Profile likelihood analysis for \code{logistiX}}
\description{
Computes the profile (penalized) likelihood at a sequence of beta values}
\usage{
\method{profile}{logistiX}(fitted, parm = 1, firth = FALSE, type = "loglike", steps = 101, betaseq = NULL, normalize = FALSE, ...)
}
\arguments{
  \item{fitted}{
A fitted \code{logistiX} object}
  \item{parm}{
Variable number for which the profile likelihood should be computed}
  \item{firth}{
if \code{TRUE}, requests penalized profile likelihood analysis}
  \item{type}{
one of \code{c("loglike","like"}. By default, returns profile log likelihood; otherwise, profile likelihood}
  \item{steps}{
  number of steps for profiling (default=101).
}
  \item{betaseq}{
  Optional: specifiy the steps at which the profile likelihood should be computed. If left \code{NULL}, the program will 
  define a range on its own slightly extending the estimated confidence limits.
}
  \item{normalize}{
  If \code{TRUE}, will normalize the profile to be 0 (for log likelihood) or 1 (for likelihood) at the maximum}
  \item{...}{
  additional argment(s) for methods.
  }
}

\details{
The profile likelihod can be checked for symmetry to see if normal approximation would be appropriate, or to evaluate whether the likelihood is flat around the estimated value.
The beta values used for computation of the profile likelihood are generated automatically, but could also be user-supplied.
}
\value{
A \code{steps} x 2 matrix of beta and profile likelihood values
}
\references{
Heinze G and Puhr R (2010). Bias-reduced and separation-proof conditional logistic regression with small or sparse data sets. Statistics in Medicine 29:770-777
}
\author{
Georg Heinze and Tobias Ladner
}
\note{
Not implemented in LogXact.
}


\seealso{
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

  # overlay unpenalized and penalized log likelihood of variable 1:
  plot(profile(test1, 1, normalize=TRUE), type="l")
  lines(profile(test1, 1, firth=TRUE, normalize=TRUE), lty=2)
  legend("topright",lty=c(1,2),legend=c("Unpenalized", "Penalized"))
}
\keyword{regression}
\keyword{models}
