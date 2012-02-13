\name{maize.model2}
\alias{maize.model2}
\title{The basic Maize model function with a litle change to use with maize.simul}
\usage{
  maize.model2(param, weather, sdate, ldate)
}
\arguments{
  \item{param}{: a vector of parameters}

  \item{weather}{: weather data.frame for one single year}

  \item{sdate}{: sowing date}

  \item{ldate}{: last date}
}
\value{
  data.frame with daily TT, LAI,B
}
\description{
  The basic Maize model function with a litle change to use
  with maize.simul
}

