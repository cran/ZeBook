#' @name Bean
#' @title Bean gene-based models dataset
#' @description Genetic data for the common bean (Phaseolus vulgaris L) that was based on a population
#' created by C. E. Vallejos (personal communication; also see Bhakta et al., 2015, 2017)
#' by crossing two widely-differing cultivars of bean (Calima, an Andean type, with Jamapa,
#' a Mesoamerican type). 
#' Bean$marker : Bean Marker Data
#' Bean$MET : Weather data on 5 Locations
#' Bean$QTL : QTL Data
#' Bean$modelpar : Dynamic Model parameters
#' @docType data
#' @usage Bean
#' @format a \code{List} including 4 data.frame  Bean$marker, Bean$MET, Bean$QTL, Bean$modelpar.
#' @source C. E. Vallejos (personal communication) and Bhakta et al. (2015, 2017).
#' Bhakta, M. S., V. A. Jones, C. E. Vallejos. 2015. Punctuated distribution of recombination hotspots and demarcation of pericentromeric regions in Phaseolus vulgaris L.  PLoS ONE 10(1): https://doi.org/10.1371/journal.pone.0116822 
#' Bhakta, M. S., S. A. Gezan, J. A. Clavijo Michelangeli, M. Carvalho, L. Zhang, J. W. Jones, K. J. Boote, M. J. Correll, J. Beaver, J. M. Osorno, R. Colbert, I. Rao, S. Beebe, A. Gonzalez, J. Ricaurte, and C. E.do Vallejos, 2017. A predictive model for time-to-flower in the common bean based on  QTL and environmental variables. G3: Genes, Genomes, Genetics 7(12) 3901-3912. https://doi.org/10.1534/g3.117.300229.
#' @examples 
#' # show the maker of JC1 to JC9 values for both parents (JAM and CAL)
#' # and 5 cRILS (RIJC001 to RIJC005)
#' Bean$marker[2:8,1:10]
#' # show the first value of weather data
#' head(Bean$MET)
#' # show the value of QTL
#' Bean$QTL[4:10,1:10]
#' # show the value of
#' Bean$modelpar
NULL
