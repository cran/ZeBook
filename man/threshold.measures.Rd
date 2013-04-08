\name{threshold.measures}
\alias{threshold.measures}
\title{TO COMPLETE Computation of threshold.measures}
\usage{
  threshold.measures(Yobs, Ypred, p, d, units = "")
}
\arguments{
  \item{Ypred}{: prediction values from the model}

  \item{Yobs}{: observed values}

  \item{p}{: TO COMPLETE}

  \item{d}{: TO COMPLETE}

  \item{units}{: units}
}
\value{
  data.frame with the different evaluation criteria
}
\description{
  TO COMPLETE Computation of threshold.measures
}
\examples{
# observed and simulated values
obs<-c(78,110,92,75,110,108,113,155,150)
sim<-c(126,126,126,105,105,105,147,147,147)
threshold.measures(obs,sim,80,1.0)
}

