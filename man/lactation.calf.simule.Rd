\name{lactation.calf.simule}
\alias{lactation.calf.simule}
\title{Wrapper function to run the Lactation model for multiple sets of parameter values}
\usage{
  lactation.calf.simule(X, duration, dt)
}
\arguments{
  \item{X}{: parameter matrix}

  \item{duration}{: duration of simulation}

  \item{dt}{: time step}
}
\value{
  data.frame with : number of paramter vector (line number
  from X), week, CS, M, Mmoy, RM, day, week
}
\description{
  Wrapper function to run the Lactation model for multiple
  sets of parameter values
}

