\name{watbal.update}
\alias{watbal.update}
\title{Water Balance model function - Update}
\usage{
  watbal.update(WAT0, RAIN, ETr, param, WP, FC)
}
\arguments{
  \item{WAT0}{: Water at the beginning of the day (mm).}

  \item{RAIN}{: Rainfall of day (mm)}

  \item{ETr}{: Evapotranspiration of day (mm)}

  \item{param}{: a vector of parameters}

  \item{FC}{: Water content at field capacity (cm^3.cm^-3)}

  \item{WP}{: Water content at wilting Point (cm^3.cm^-3)}
}
\value{
  WAT1 : Water at the beginning of the day+1 (mm).
}
\description{
  Water Balance model function - Update
}

