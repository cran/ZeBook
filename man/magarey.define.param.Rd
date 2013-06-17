\name{magarey.define.param}
\alias{magarey.define.param}
\title{Define values of the parameters for the Magarey model}
\usage{
  magarey.define.param(species = "unkown")
}
\arguments{
  \item{species}{: name of a species. By default, value for
  an "unkown" species are given. Other possibility are
  "G.citricarpa" or "Kaki.fungus"}
}
\value{
  matrix with parameter values (nominal, binf, bsup)
}
\description{
  Define values of the parameters for the Magarey model
}
\examples{
magarey.define.param(species="G.citricarpa")
magarey.define.param(species="Kaki.fungus")
}

