\name{maize_cir_rue.model}
\alias{maize_cir_rue.model}
\title{The Maize model function with rue (and CumInt)}
\usage{
  maize_cir_rue.model(Tbase, RUE_max, K, alpha, LAImax,
    TTM, TTL, weather, sdate, ldate)
}
\arguments{
  \item{Tbase,RUE_max,K,alpha,LAImax,TTM,TTL}{: 7
  parameters, weather table, sowing date and last date}

  \item{weather}{: weather data.frame for one single year}

  \item{sdate}{: sowing date}

  \item{ldate}{: last date}
}
\value{
  data.frame with daily TT, LAI,B
}
\description{
  The Maize model function with rue (and CumInt)
}

