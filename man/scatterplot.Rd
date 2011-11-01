\name{scatterplot}
\alias{scatterplot}
\title{Create scatter plot with imputed survival times}
\description{
Create scatter plot with imputed survival times.
}
\usage{scatterplot(formula,data,horizon,plot=TRUE,xlab)}
\arguments{
  \item{formula}{Formula for prediction model to be used as in
  \code{\link[survival:coxph]{coxph}}}
  \item{data}{Data set in which to interpret the formula}
  \item{horizon}{The horizon, maximum value to be imputed in case of censored
  observations; default is 1.05 times largest event time}
  \item{plot}{Should the tolerance plot actually be plotted? Default is
  \code{TRUE}}
  \item{xlab}{Label for x-axis}
}
\details{
Imputation is used for censored survival times.}
\value{
A data frame with columns
  \item{x}{Predictor (centered at zero)}
  \item{imputed}{(Imputed) survival time}
and with attribute \code{"horizon"} (copied from input or default).}
\references{
Royston P (2001), The lognormal distribution as a model for survival time
in cancer, with an emphasis on prognostic factors, Statistica Neerlandica
55, 89-104.

van Houwelingen HC, Putter H (2011). Dynamic Prediction in Clinical Survival Analysis.
Chapman & Hall.
}
\author{Hein Putter \email{H.Putter@lumc.nl}}
\examples{
data(ova)
scatterplot(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova)
}
\keyword{univar}
