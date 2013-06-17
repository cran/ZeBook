\name{exponential.model.ie}
\alias{exponential.model.ie}
\title{The Exponential growth model of dynamic of population - with improved Euler integration}
\usage{
  exponential.model.ie(a, Y0, duration = 40, dt = 1)
}
\arguments{
  \item{a}{: growth rate}

  \item{Y0}{: initial condition}

  \item{duration}{: duration of simulation}

  \item{dt}{: time step for integration}
}
\value{
  data.frame with Y for each time step
}
\description{
  The Exponential growth model of dynamic of population -
  with improved Euler integration
}
\seealso{
  \code{\link{verhulst.update}} for the update function of
  the Verhulst model.
}

