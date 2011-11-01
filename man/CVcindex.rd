\name{CVcindex}
\alias{CVcindex}
\title{Calculate cross-validated c-index}
\description{
This function calculates cross-validated versions of Harrell's c-index.}
\usage{CVcindex(formula,data,type="single",matrix=FALSE)}
\arguments{
  \item{formula}{Formula for prediction model to be used as in
  \code{\link[survival:coxph]{coxph}}}
  \item{data}{Data set in which to interpret the formula}
  \item{type}{One of \code{"single"}, \code{"pair"} or \code{"fullpairs"}. For
  \code{"single"} (default), the prognostic index Z_i is replaced by Z_{i,(-i)},
  for \code{"pair"}, two assessments of concordance are made for each pair
  (i,j), one using Z_{i,(-i)} and Z_{j,(-i)}, the other using Z_{i,(-j)} and
  Z_{j,(-j)}, for \code{"fullpairs"}, each of the possible pairs is left out
  and comparison is based on Z_{i,(-i,-j)} and Z_{j,(-i,-j)}}
  \item{matrix}{if \code{TRUE}, the matrix of cross-validated prognostic indices
  is also returned; default is \code{FALSE}}
}
\value{
A list with elements
\item{concordant}{The number of concordant pairs}
\item{total}{The total number of pairs that can be evaluated}
\item{cindex}{The cross-validated c-index}
\item{matrix}{Matrix of cross-validated prognostic indices (only if argument
  \code{matrix} is \code{TRUE}}
and with attribute \code{"type"} as given as input.
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
# Real thing takes a long time, so on a smaller data set
ova2 <- ova[1:80,]
# Actual c-index
cindex(Surv(tyears,d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova2)
# Cross-validated c-indices
CVcindex(Surv(tyears,d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova2)
CVcindex(Surv(tyears,d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova2, type="pair")
CVcindex(Surv(tyears,d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova2, type="fullpairs")
}
\keyword{univar}