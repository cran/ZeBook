\name{magarey.model}
\alias{magarey.model}
\title{The Magarey model}
\usage{
  magarey.model(T, Tmin, Topt, Tmax, Wmin, Wmax)
}
\arguments{
  \item{T}{: input variable. Either a scalar or a vector
  (for a weather serie).}

  \item{Tmin}{: parameter of minimal temperature for
  infection (degC)}

  \item{Topt}{: parameter of optimal temperature for
  infection (degC)}

  \item{Tmax}{: parameter of maximal temperature for
  infection (degC)}

  \item{Wmin}{: parameter of minimal wetness duration for
  infection (hour)}

  \item{Wmax}{: parameter of maximal wetness duration for
  infection (hour)}
}
\value{
  Wetness duration (W, hour). Either a scalar or a vector
  depending on T.
}
\description{
  Generic model of infection for foliar diseases caused by
  fungi (from Magarey et al.,2005).
}

