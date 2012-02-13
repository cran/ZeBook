\name{weed.update}
\alias{weed.update}
\title{The weed update function}
\usage{
  weed.update(d, S, SSBa, DSBa, Soil, Crop, Herb, param)
}
\arguments{
  \item{d}{: weed density at seed emergence (plants/m2) -
  value for year}

  \item{S}{: seed production per m2 - value for year}

  \item{SSBa}{: surface seedbank after tillage (grains/m2)
  - value for year}

  \item{DSBa}{: deep seedbank after tillage (grains/m2) -
  value for year}

  \item{Soil}{: value for soil decision (1 or 0)}

  \item{Crop}{: value for crop decision (1 or 0)}

  \item{Herb}{: value for herbicide traitment decision (1
  or 0)}

  \item{param}{: vector of the 16 parameters}
}
\value{
  a vector with values of state variables for year+1
}
\description{
  TOCOMPLETE
}

