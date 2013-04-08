\name{magarey.simule}
\alias{magarey.simule}
\title{Wrapping function to run an virtual experimental design for Magarey model}
\usage{
  magarey.simule(X, T, all = FALSE)
}
\arguments{
  \item{X}{: parameter matrix}

  \item{T}{: input variable, temperature}

  \item{all}{: if you want a matrix combining X and output}
}
\value{
  a table with wetness duration (W) for each parameter
  vector
}
\description{
  Example magarey.simule(magarey.define.param(),15)
}

