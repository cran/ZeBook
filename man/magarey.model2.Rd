\name{magarey.model2}
\alias{magarey.model2}
\title{The Magarey model, taking in argument a vector of param}
\usage{
  magarey.model2(T, param)
}
\arguments{
  \item{param}{: parameters}

  \item{T}{: input variable. Either a scalar or a vector
  (for a weather series).}
}
\value{
  W : Wetness duration (hour). Either a scalar or a vector
  depending on T.
}
\description{
  Generic model of infection for foliar diseases caused by
  fungi (from Magarey et al.,2005).
}

