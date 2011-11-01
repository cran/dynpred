\name{toleranceplot}
\alias{toleranceplot}
\title{Create a tolerance plot}
\description{
Create a tolerance plot according to 
}
\usage{toleranceplot(formula,data,coverage=0.8,horizon,plot=TRUE,xlab)}
\arguments{
  \item{formula}{Formula for prediction model to be used as in
  \code{\link[survival:coxph]{coxph}}}
  \item{data}{Data set in which to interpret the formula}
  \item{coverage}{The coverage for the tolerance intervals (default is 0.8)}
  \item{horizon}{The horizon, maximum value to be imputed in case of censored
  observations; default is 1.05 times largest event time}
  \item{plot}{Should the tolerance plot actually be plotted? Default is
  \code{TRUE}}
  \item{xlab}{Label for x-axis}
}
\details{
Warnings will be issued each time the survival curve corresponding to a value of
x never goes below (1-coverage)/2; these warnings may be ignored. 
}
\value{
A data frame with columns
  \item{x}{Predictor (centered at zero)}
  \item{lower}{Lower bound of tolerance interval}
  \item{upper}{Upper bound of tolerance interval}
and with attributes \code{"coverage"} and \code{"horizon"} (copied from input or
default).}
\references{
Henderson R, Jones M & Stare J (2001), Accuracy of point predictions in survival
analysis, Statistics in Medicine 20, 3083-3096.

van Houwelingen HC, Putter H (2011). Dynamic Prediction in Clinical Survival Analysis.
Chapman & Hall.
}
\author{Hein Putter \email{H.Putter@lumc.nl}}
\examples{
data(ova)
toleranceplot(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova)
}
\keyword{univar}
