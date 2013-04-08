\name{epirice.multi.simul}
\alias{epirice.multi.simul}
\title{Wrapping function to run an virtual experimental design for EPIRICE model}
\usage{
  epirice.multi.simul(param, multi.simul, all = FALSE)
}
\arguments{
  \item{param}{: a vector of parameters}

  \item{multi.simul}{: matrix of n row definition of input
  variable : site, year and date of transplantation.}

  \item{all}{: if you want a matrix combining multi.simul
  and output (default = FALSE)}
}
\value{
  matrix with AUDPC for each input vector
}
\description{
  Wrapping function to run an virtual experimental design
  for EPIRICE model
}
\seealso{
  \code{\link{epirice.model}}
}

