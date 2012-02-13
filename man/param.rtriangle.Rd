\name{param.rtriangle}
\alias{param.rtriangle}
\title{Building a random plan with a triangle distribution}
\usage{
  param.rtriangle(model.factors, N)
}
\arguments{
  \item{model.factors}{: matrix defining nominal, minimal
  (binf), maximal values (bsup) for a set of p parameters}

  \item{N}{: size of sample}
}
\value{
  parameter matrix of dim = (N, p)
}
\description{
  according to nominal, minimal and maximal values defined
  in a model.factors matrix
}

