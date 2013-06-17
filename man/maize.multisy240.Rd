\name{maize.multisy240}
\alias{maize.multisy240}
\title{Wrapper function to run Maize model for multiple sets of input variables (site-year) and give Biomass at day240.}
\usage{
  maize.multisy240(param, liste_sy, sdate, ldate)
}
\arguments{
  \item{param}{: a vector of parameters}

  \item{liste_sy}{: vector of site-year}

  \item{sdate}{: sowing date}

  \item{ldate}{: last date}
}
\value{
  mean biomass at day=240
}
\description{
  Wrapper function to run Maize model for multiple sets of
  input variables (site-year) and give Biomass at day240.
}
\examples{
maize.multisy240(maize.define.param()["nominal",],c("18-2006","64-2004") , sdate=100, ldate=250)
}

