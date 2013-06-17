\name{carcass.EMI.multi}
\alias{carcass.EMI.multi}
\title{Wrapper function to run Carcass model on several animals with different conditions}
\usage{
  carcass.EMI.multi(param, list_individuals, energy,
    init_condition)
}
\arguments{
  \item{param}{: a vector of parameters}

  \item{list_individuals}{: list of individuals}

  \item{energy}{: Metabolizable energy available for all
  individuals}

  \item{init_condition}{: initial condition for all
  individuals}
}
\value{
  data.frame with id, ration ,duration, day,
  PV,ProtC,ProtNC,LipC,LipNC
}
\description{
  Wrapper function to run Carcass model on several animals
  with different conditions
}

