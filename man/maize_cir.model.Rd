\name{maize_cir.model}
\alias{maize_cir.model}
\title{the Maize model function with additionnal state varaiable CumInt}
\usage{
  maize_cir.model(Tbase, RUE, K, alpha, LAImax, TTM, TTL,
    weather, sdate, ldate)
}
\arguments{
  \item{Tbase,RUE,K,alpha,LAImax,TTM,TTL}{: 7 parameters}

  \item{weather}{: weather data.frame for one single year}

  \item{sdate}{: sowing date}

  \item{ldate}{: last date}
}
\value{
  data.frame with daily TT, LAI,B
}
\description{
  the Maize model function with additionnal state varaiable
  CumInt
}

