\name{cindex}
\alias{cindex}
\title{Calculate Harrell's c-index}
\description{
This function calculates Harrell's c-index.}
\usage{cindex(formula,data)}
\arguments{
  \item{formula}{Formula for prediction model to be used as in
  \code{\link[survival:coxph]{coxph}}}
  \item{data}{Data set in which to interpret the formula}
}
\value{
A list with elements
\item{concordant}{The number of concordant pairs}
\item{total}{The total number of pairs that can be evaluated}
\item{cindex}{Harrell's c-index}
}
\references{
Harrell FE, Lee KL & Mark DB (1996), Multivariable prognostic models: issues in
developing models, evaluating assumptions and adequacy, and measuring and reducing
errors, Statistics in Medicine 15, 361-387.

van Houwelingen HC, Putter H (2011). Dynamic Prediction in Clinical Survival Analysis.
Chapman & Hall.
}
\author{Hein Putter \email{H.Putter@lumc.nl}}
\examples{
data(ova)
cindex(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova)
}
\keyword{univar}