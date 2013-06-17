\name{lactation.model.machine2}
\alias{lactation.model.machine2}
\title{The Lactation model for use with lactation.machine.simule}
\usage{
  lactation.model.machine2(param, duration, dt, CSi, Mi)
}
\arguments{
  \item{param}{: a vector of parameters containning
  (cu,kdiv,kdl,kdh,km,ksl,kr,ks,ksm,mh,mm,p,mum,rma,t1,t2,t3,t4,t5,t6)(see
  lactation.model.machine)}

  \item{duration}{: duration of simulation}

  \item{dt}{: time step}

  \item{CSi}{: initial Number of secretory cells}

  \item{Mi}{: initial Quantity of milk in animal (kg)}
}
\value{
  data.frame with CS, M, Mmoy, RM, day, week
}
\description{
  see lactation.calf.model for model description.
}

