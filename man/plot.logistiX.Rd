\name{plot.logistiX}
\alias{plot.logistiX}
\title{
Plot method for logistiX objects
}
\description{
Plots the exact null or alternative distribution of the sufficient statistic for a variable of interest. Computes the likelihood or penalized likelihood
for the observed value of that variable.
}
\usage{
\method{plot}{logistiX}(x, var, beta = 0, plot = TRUE, firth = FALSE, ...)
}
\arguments{
  \item{x}{
  a \code{logistiX} object
  }
  \item{var}{
  The variable number, for which the distribution of sufficient statistic should be plotted.
  }
  \item{beta}{
  The value of the regression coefficient for which the exact distribution should be computed. (default=0, will plot the null distribution)
}
  \item{plot}{
  Set to \code{TRUE} for plotting (default=\code{TRUE}). Set to \code{FALSE} to suppress plotting.
 }
  \item{firth}{
  Set to \code{TRUE} for computing the penalized instead of the ordinary likelihood.}
  \item{...}{further arguments}
}
\details{
This method can be used to visualize the exact distribution of the sufficient statistic of a variable of interest, and/or it can be used to compute the likleihood
or penalized likelihood for the sufficient statistic of that variable.

The method can also display the distribution under a non-0 beta coefficient.
}
\value{
A list with entries:
\item{\code{dist}}{ a data frame with columns 
  \code{t.stat}: the value of the sufficient statistic
  \code{score}: the probability for the sufficient statistic)
  \code{penalized.score}: the penalized likelihood of each value of the sufficient statistic
  }
\item{\code{likelihood}}{the value of the (penalized) likelihood at the observed value of the sufficient statistic}
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

  plot(test1,1) # plots the distribution for the first variable
  plot(test1,1, plot=FALSE, firth=TRUE)$likelihood # returns the penalized null likelihood for the first variable
}  

\keyword{regression}
\keyword{models}
