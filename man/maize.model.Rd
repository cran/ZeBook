\name{maize.model}
\alias{maize.model}
\title{The basic Maize model function}
\usage{
  maize.model(Tbase, RUE, K, alpha, LAImax, TTM, TTL,
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
  \strong{Model description.} This model is a dynamic model
  of crop growth for Maize cultivated in potential
  conditions. The crop growth is represented by three state
  variables, leaf area per unit ground area (leaf area
  index, LAI), total biomass (B) and cumulative thermal
  time since plant emergence (TT). It is based on key
  concepts included in most crop models, at least for the
  "potential production" part. In fact, this model does not
  take into account any effects of soil water, nutrients,
  pests, or diseases,... The tree state variables are
  dynamic variables depending on days after emergence:
  TT(day), B(day), and LAI(day). The model has a time step
  dt of one day.\cr The model is defined by a few
  equations, with a total of seven parameters for the
  described process. \cr (1) \eqn{TT(day + 1) = TT(day ) +
  dTT(day)}
}

