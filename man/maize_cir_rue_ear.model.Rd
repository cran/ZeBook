\name{maize_cir_rue_ear.model}
\alias{maize_cir_rue_ear.model}
\title{The Maize model function with CumInt, RUE and ear growth}
\usage{
  maize_cir_rue_ear.model(Tbase, RUE_max, K, alpha, LAImax,
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
  The Maize model function with CumInt, RUE and ear growth
}

