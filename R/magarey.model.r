################################################################################
# "Working with dynamic models for agriculture"
# R script for pratical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2010-08-09
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
################################ FUNCTIONS #####################################
#' @title The Magarey model
#' @description Generic model of infection for foliar diseases caused by fungi (from Magarey et al.,2005).
#' @param T : input variable. Either a scalar or a vector (for a weather serie).
#' @param Tmin : parameter of minimal temperature for infection (degC)
#' @param Topt : parameter of optimal temperature for infection (degC)
#' @param Tmax : parameter of maximal temperature for infection (degC)
#' @param Wmin : parameter of minimal wetness duration for infection (hour)
#' @param Wmax : parameter of maximal wetness duration for infection (hour)
#' @return Wetness duration (W, hour). Either a scalar or a vector depending on T. 
#' @export
magarey.model<-function(T, Tmin, Topt, Tmax, Wmin, Wmax){
    fT<-((Tmax-T)/(Tmax-Topt))*(((T-Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt)))
	W <- Wmin/fT
	W[W>Wmax | T<Tmin | T>Tmax]<-Wmax 
	return(W)
}
#' @title The Magarey model, taking in argument a vector of param
#' @description Generic model of infection for foliar diseases caused by fungi (from Magarey et al.,2005).
#' @param param : parameters 
#' @param T : input variable. Either a scalar ou a vector (for a weather serie). 
#' @return W : Wetness duration (hour). Either a scalar or a vector depending on T.  
#' @export
magarey.model2<-function(T, param){
    return(magarey.model(T, param["Tmin"], param["Topt"], param["Tmax"], param["Wmin"], param["Wmax"]))
}


################################################################################
#' @title Define parameter values of the model function
#' @return matrix with parameter values (nominal, binf, bsup)
#' @export
magarey.define.param <- function()
{
# nominal, binf, bsup
# Tmin  : the minimal temperature for infection (°C)
Tmin = c(NA, 0.8, 1.2)
# Topt  : the optimal temperature for infection (°C)
Topt = c(NA, 20, 30)
# Tmax  : the maximal temperature for infection (°C)
Tmax = c(NA, 24, 36)
# Wmin : minimal wetness duration for infection (hour)
Wmin = c(NA, 38.4, 57.6)
# Wmax : maximal wetness duration for infection (hour)
Wmax = c(NA, 115.2, 172.8)

param<-data.frame(Tmin,Topt,Tmax, Wmin, Wmax)
row.names(param)<-c("nominal","binf","bsup")
return(as.matrix(param))
}
################################################################################
#' @title Wrapping fonction to run an virtual experimental design for Magarey model
#' @param X : parameter matrix
#' @param T : input variable, temperature
#' @param all : if you want a matrix combining X and output
#' @return a table with wetness duration (W) for each parameter vector
#' @export
#' @description Example magarey.simule(magarey.define.param(),15)
magarey.simule <- function(X, T,  all=FALSE){
Y <- apply(X,1,function(param) magarey.model2(T, param))
if(all) Y = cbind(X,W = Y)
return(as.matrix(Y))
}


################################################################################
# End of file