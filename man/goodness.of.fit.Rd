\name{goodness.of.fit}
\alias{goodness.of.fit}
\title{Computation of classical evaluation criteria function - NEW}
\usage{
  goodness.of.fit(Yobs, Ypred, draw.plot = FALSE)
}
\arguments{
  \item{Ypred}{: prediction values from the model}

  \item{Yobs}{: observed values}

  \item{draw.plot}{: draw evaluation plot}
}
\value{
  data.frame with the different evaluation criteria
}
\description{
  Computation of classical evaluation criteria function -
  NEW
}
\examples{
# observed and simulated values
obs<-c(78,110,92,75,110,108,113,155,150)
sim<-c(126,126,126,105,105,105,147,147,147)
goodness.of.fit(obs,sim,draw.plot=TRUE)
}

