#' @name weather_EuropeEU
#' @title Weather Series for Europe (EU) from NASA POWER Agroclimatology
#'
#' @description
#' A contemporary daily climate dataset for Europe, covering the period from
#' 1st January 2001 to 31 December 2010 (10 complete years of data).
#' This dataset focuses on regions of Europe with an elevation less than 500 meters.
#'
#' The dataset was extracted from the NASA Langley Research Center POWER Project,
#' which provides agroclimatology datasets (Chandler et al., 2004).
#' It was funded through the NASA Earth Science Directorate Applied Science Program.
#'
#' This climate dataset contains daily estimates of:
#' \itemize{
#'   \item Precipitation
#'   \item Mean, minimum, and maximum temperature
#'   \item Relative humidity
#'   \item Dew point
#'   \item Solar radiation
#'   \item Wind speed
#' }
#'
#' The data has global coverage at one-degree resolution (approximately 111 km at the equator).
#'
#' The NASA POWER agroclimatology data are derived from multiple sources:
#' \itemize{
#'   \item Solar radiation: satellite observations
#'   \item Meteorological data: Goddard Earth Observing System global assimilation model version 4 (GEOS-4)
#'   \item Precipitation: Global Precipitation Climate Project and Tropical Rainfall Measurement Mission
#' }
#'
#' A full description is available at \url{https://power.larc.nasa.gov}.
#'
#' Elevation (Altitude) was retrieved from the Aster Global Digital Elevation Model
#' using the webservice api.geonames.org/astergdem.
#' Samples are approximately 30m x 30m, covering latitudes between 83 degrees N and 65 degrees S.
#' Results provide a single number for elevation in meters (Aster GDEM).
#' Ocean areas are masked as "no data" and assigned a value of -9999.
#'
#' @docType data
#' @usage data(weather_EuropeEU)
#'
#' @format
#' A \code{RangedData} instance with 1 row per day and the following variables:
#' \describe{
#'   \item{SRAD}{Daily Insolation Incident On A Horizontal Surface (MJ/m2/day)}
#'   \item{T2M}{Average Air Temperature At 2 m Above The Surface Of The Earth (degrees C)}
#'   \item{TMIN}{Minimum Air Temperature At 2 m Above The Surface Of The Earth (degrees C)}
#'   \item{TMAX}{Maximum Air Temperature At 2 m Above The Surface Of The Earth (degrees C)}
#'   \item{RH2M}{Relative Humidity At 2 m (percent)}
#'   \item{TDEW}{Dew/Frost Point Temperature At 2 m (degrees C)}
#'   \item{RAIN}{Average Precipitation (mm/day)}
#'   \item{WIND}{Wind Speed At 10 m Above The Surface Of The Earth (m/s)}
#' }
#'
#' @source
#' \url{https://power.larc.nasa.gov/}
#' \url{https://asterweb.jpl.nasa.gov/gdem.asp}
#' \url{https://www.geonames.org/about.html}
#'
#' @examples
#' data(weather_EuropeEU)
#' head(weather_EuropeEU)
NULL