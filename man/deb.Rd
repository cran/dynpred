\name{deb}
\alias{deb}
\title{Debugging function}
\description{
A simple but useful debugging function. It first announces the object to printed and then prints it.}
\usage{deb(x, method=c("print","cat"))}
\arguments{
  \item{x}{The object to be printed}
  \item{method}{The method for printing \code{x}. Default is \code{"print"}, which
  uses \code{\link{print}} for printing; \code{"cat"} uses \code{\link{cat}} for
  printing. The latter is useful for short objects (scalar and vectors), the former
  for more structured objects (data frames, matrices, lists etc).}
}
\author{Hein Putter \email{H.Putter@lumc.nl}}
\examples{
tm <- c(0.2,0.5,1,1.2,1.8,4)
ta <- 2*tm
dfr <- data.frame(time=tm, stepf=ta)
deb(dfr, method="print")
deb(tm, method="cat")
}
\keyword{univar}
