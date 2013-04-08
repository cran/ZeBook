\name{maize.simule_multisy240}
\alias{maize.simule_multisy240}
\title{Wrapping function to run maize model and compute mean Biomass at day240 for a list of site-year, for a matrix of parameter (design).}
\usage{
  maize.simule_multisy240(X, liste_sy, sdate, ldate,
    all = FALSE)
}
\arguments{
  \item{X}{: matrix of n row vectors of 7 parameters}

  \item{liste_sy}{: vector of site-year}

  \item{sdate}{: sowing date}

  \item{ldate}{: last date}

  \item{all}{: if you want a matrix combining X and output
  (default = FALSE)}
}
\value{
  a matrix of mean biomass at day=240 for all combinations
  of parameters of X
}
\description{
  Wrapping function to run maize model and compute mean
  Biomass at day240 for a list of site-year, for a matrix
  of parameter (design).
}
\examples{
maize.simule_multisy240(maize.define.param(),c("18-2006","64-2004") , sdate=100, ldate=250, all=FALSE)
}

