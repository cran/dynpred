\name{AUC}
\alias{AUC}
\title{Calculate AUC(t) curve}
\description{
Calculate model-free curve of Area Under the Curve values over time, based on the
dynamic/xxxx AUC of Heagerty and Zheng.}
\usage{AUC(formula,data,plot=TRUE)}
\arguments{
  \item{formula}{Formula for prediction model to be used as in
  \code{\link[survival:coxph]{coxph}}}
  \item{data}{Data set in which to interpret the formula}
  \item{plot}{Determines whether the AUC function should be plotted (if \code{TRUE} (default))
  along with a \code{\link{lowess}} curve or not (if \code{FALSE})}
}
\value{
A list with elements
\item{AUCt}{A data frame with time t in column \code{time} and AUC(t) in
column \code{AUC}}
\item{AUC}{The AUC(t) weighted by Y(t)-1, with Y(t) the number at risk at t; this coincides
with Harrell's c-index}
}
\references{
Harrell FE, Lee KL & Mark DB (1996), Multivariable prognostic models: issues in
developing models, evaluating assumptions and adequacy, and measuring and reducing
errors, Statistics in Medicine 15, 361-387.

Heagerty PJ & Zheng Y (2005), Survival model predictive accuracy and ROC curves,
Biometrics 61, 92-105.

van Houwelingen HC, Putter H (2011). Dynamic Prediction in Clinical Survival Analysis.
Chapman & Hall.
}
\author{Hein Putter \email{H.Putter@lumc.nl}}
\examples{
data(ova)
AUC(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova)
}
\keyword{univar}
