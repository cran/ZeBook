\name{verhulst.model}
\alias{verhulst.model}
\title{the Verhulst model function - integration loop}
\usage{
  verhulst.model(a, k, Y0, duration)
}
\arguments{
  \item{a}{: growth rate}

  \item{k}{: capacity}

  \item{Y0}{: initial condition}

  \item{duration}{: duration of simulation}
}
\value{
  data.frame with daily Y
}
\description{
  the Verhulst model function - integration loop
}
\seealso{
  \code{\link{verhulst.update}} for the update function of
  the Verhulst model.
}

