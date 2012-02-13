################################################################################
# "Working with dynamic models for agriculture"
# R script for pratical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2010-08-09
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
################################ FUNCTIONS #####################################
#' @title The basic Maize model function
#' @description \strong{Model description.}
#' This model is a dynamic model of crop growth for Maize cultivated in potential conditions.
#' The crop growth is represented by three state variables, leaf area per unit ground area (leaf area index, LAI), total biomass (B) and cumulative thermal time since plant emergence (TT). It is based on key concepts included in most crop models, at least for the "potential production" part. In fact, this model does not take into account any effects of soil water, nutrients, pests, or diseases,... 
#' The tree state variables are dynamic variables depending on days after emergence: TT(day), B(day), and LAI(day). The model has a time step dt of one day.\cr
#' The model is defined by a few equations, with a total of seven parameters for the described process.
#' \cr (1) \eqn{TT(day  +  1)  =  TT(day )  +  dTT(day)}
#' @param Tbase,RUE,K,alpha,LAImax,TTM,TTL : 7 parameters
#' @param weather : weather data.frame for one single year
#' @param sdate : sowing date
#' @param ldate : last date
#' @return data.frame with daily TT, LAI,B
#' @export
maize.model<-function(Tbase,RUE,K,alpha,LAImax,TTM,TTL,weather,sdate,ldate)
    {
    # Initialize variables
    # 3 states variables, as 3 vectors initialized to NA
    # TT : temperature sum (°C.d)
    TT<-rep(NA,ldate)
    # B : Biomass (g/m2)
    B<-rep(NA,ldate)
    # LAI : Leaf Area Index (m2 leaf/m2 soil)
    LAI<-rep(NA,ldate)
    
    # Initialize state variables when sowing on day "sdate"   
    TT[sdate]<- 0
    B[sdate]<- 1
    LAI[sdate]<- 0.01
    
    # Simulation loop
    for (day in sdate:(ldate-1))
        {
        # Calculate rates of change of state variables (dTT, dB, dLAI)
        dTT <- max((weather$Tmin[day]+weather$Tmax[day])/2-Tbase, 0) 
        if (TT[day]<=TTM) {dB <- RUE*(1-exp(-K*LAI[day]))*weather$I[day]}
        else {dB <- 0}
        if (TT[day]<=TTL) {dLAI <- alpha*dTT*LAI[day]*max(LAImax-LAI[day],0)}
        else {dLAI <-0 }
        
        # Uptade state variables 
        TT[day+1]<- TT[day] + dTT 
        B[day+1]<- B[day] + dB
        LAI[day+1]<- LAI[day] + dLAI   
        }
        # End simulation loop
    return(data.frame(day=sdate:ldate,TT=TT[sdate:ldate],LAI=LAI[sdate:ldate],B=B[sdate:ldate]))    
    }
################################################################################
#' @title The basic Maize model function with a litle change to use with maize.simul
#' @param param : a vector of parameters
#' @param weather : weather data.frame for one single year
#' @param sdate : sowing date
#' @param ldate : last date
#' @return data.frame with daily TT, LAI,B
#' @export
maize.model2<-function(param, weather,sdate,ldate)
{
# 7 Parameter values of the model, read from the param vector
Tbase <- param["Tbase"]
RUE <- param["RUE"]
K <- param["K"]
alpha <- param["alpha"]
LAImax <- param["LAImax"]
TTM <- param["TTM"]
TTL <- param["TTL"]
    # Initialize variables
    # 3 states variables, as 3 vectors initialized to NA
    # TT : temperature sum (°C.d)
    TT<-rep(NA,ldate)
    # B : Biomass (g/m2)
    B<-rep(NA,ldate)
    # LAI : Leaf Area Index (m2 leaf/m2 soil)
    LAI<-rep(NA,ldate)

    # Initialize state variables when sowing on day "sdate"
    TT[sdate]<- 0
    B[sdate]<- 1
    LAI[sdate]<- 0.01

    # Simulation loop
    for (day in sdate:(ldate-1))
        {
        # Calculate rates of change of state variables (dTT, dB, dLAI)
        dTT <- max((weather$Tmin[day]+weather$Tmax[day])/2-Tbase, 0)
        if (TT[day]<=TTM) {dB <- RUE*(1-exp(-K*LAI[day]))*weather$I[day]}
        else {dB <- 0}
        if (TT[day]<=TTL) {dLAI <- alpha*dTT*LAI[day]*max(LAImax-LAI[day],0)}
        else {dLAI <-0 }

        # Uptade state variables
        TT[day+1]<- TT[day] + dTT
        B[day+1]<- B[day] + dB
        LAI[day+1]<- LAI[day] + dLAI
        }
        # End simulation loop
    return(data.frame(day=sdate:ldate,TT=TT[sdate:ldate],LAI=LAI[sdate:ldate],B=B[sdate:ldate]))
}
################################################################################
#' @title Define parameter values of the maize model
#' @return matrix with parameter values (nominal, binf, bsup)
#' @export
maize.define.param <- function()
{
# nominal, binf, bsup
# Tbase  : the baseline temperature for growth (°C)
Tbase <- c(7, 6, 8)
# RUE : radiation use efficiency (g.MJ-1)
RUE <- c(1.85,1.5,2.5)
# K : extinction coefficient (-)
K <- c(0.7,0.6,0.8)
#alpha : the relative rate of leaf area index increase for small values of leaf area index ((°C.day)-1)
alpha <- c(0.00243,0.002,0.003)
#LAImax : maximum leaf area index (m2 leaf/m2 soil)
LAImax <- c(7.0,6.0,8.0)
#TTM :  temperature sum for crop maturity (°C.day)
TTM <- c(1200,1100,1400)
#TTL : temperature sum at the end of leaf area increase (°C.day)
TTL <- c(700,600,850)
param<-data.frame(Tbase,RUE,K, alpha, LAImax, TTM, TTL)
row.names(param)<-c("nominal","binf","bsup")
return(as.matrix(param))
}
################################################################################
#' @title Wrapping fonction to run an virtual experimental design for maize model
#' @param X : matrix of n row vectors of 7 parameters
#' @param weather : weather data.frame for one single year
#' @param sdate : sowing date
#' @param ldate : last date
#' @param all : if you want a matrix combining X and output (delfaut = FALSE)
#' @return matrix with maximum biomass for each parameter vector
#' @export
maize.simule <- function(X, weather, sdate, ldate, all=FALSE){
# output : maximum biomass only
Y <- apply(X,1,function(v) max(maize.model2(v[1:7],weather, sdate, ldate)$B, na.rm=TRUE))
if(all) Y = cbind(X,B = Y)
return(as.matrix(Y))
}
################################################################################
#' @title the Maize model function with additionnal state varaiable CumInt
#' @param Tbase,RUE,K,alpha,LAImax,TTM,TTL : 7 parameters
#' @param  weather : weather data.frame for one single year
#' @param sdate : sowing date
#' @param ldate : last date
#' @return data.frame with daily TT, LAI,B
#' @export
maize_cir.model<-function(Tbase,RUE,K,alpha,LAImax,TTM,TTL,weather,sdate,ldate)
    {
    # Initialize variables
    # 3 states variables, as 3 vectors initialized to NA
    # TT : temperature sum (°C.d)
    TT<-rep(NA,ldate)
    # B : Biomass (g/m2)
    B<-rep(NA,ldate)
    # LAI : Leaf Area Index (m2 leaf/m2 soil)
    LAI<-rep(NA,ldate)
    # CumInt LAI : Cumulative intercepted radiation
	CumInt<-rep(NA,ldate) 


    # Initialize state variables when sowing on day "sdate"   
    TT[sdate]<- 0
    B[sdate]<- 1
    LAI[sdate]<- 0.01
    CumInt[sdate] = 0.0

    # Simulation loop
    for (day in sdate:(ldate-1))
        {
        # Calculate rates of change of state variables (dTT, dB, dLAI)
        dTT <- max((weather$Tmin[day]+weather$Tmax[day])/2-Tbase, 0) 
        if (TT[day]<=TTM) {dB <- RUE*(1-exp(-K*LAI[day]))*weather$I[day]}
        else {dB <- 0}
        if (TT[day]<=TTL) {dLAI <- alpha*dTT*LAI[day]*max(LAImax-LAI[day],0)}
        else {dLAI <-0 }
		
        
        # Uptade state variables 
        TT[day+1]<- TT[day] + dTT 
        B[day+1]<- B[day] + dB
        LAI[day+1]<- LAI[day] + dLAI   
		CumInt[day+1] = CumInt[day] + weather$I[day]*(1 - exp(- K * LAI[day]))
        }
        # End simulation loop
    return(data.frame(day=sdate:ldate,TT=TT[sdate:ldate],LAI=LAI[sdate:ldate],B=B[sdate:ldate],CumInt=CumInt[sdate:ldate]))    
    }

###############################################################################
#' @title Effect of the temperature on RUE for Maize
#' @param T : temperature
#' @param RUE_max : maximum value for RUE
#' @param T0 : temperature parameter
#' @param T1 : temperature parameter
#' @param T2 : temperature parameter
#' @param T3 : temperature parameter
#' @return RUE value
#' @export
maize.RUEtemp <- function(T, RUE_max,T0,T1,T2,T3)
	{
	RUE = ((T>=T0)*(T<T1))* RUE_max*(T-T0)/(T1-T0) + ((T>=T1)*(T<T2))*RUE_max + ((T>=T2)*(T<T3))*RUE_max*(T3-T)/(T3-T2)
	}
###############################################################################
#' @title The Maize model function with rue (and CumInt)
#' @param Tbase,RUE_max,K,alpha,LAImax,TTM,TTL : 7 parameters, weather table, sowing date and last date
#' @param  weather : weather data.frame for one single year
#' @param sdate : sowing date
#' @param ldate : last date
#' @return data.frame with daily TT, LAI,B
#' @export
maize_cir_rue.model<-function(Tbase,RUE_max,K,alpha,LAImax,TTM,TTL,weather,sdate,ldate)
    {
    # Initialize variables
    # 3 states variables, as 3 vectors initialized to NA
    # TT : temperature sum (°C.d)
    TT<-rep(NA,ldate)
    # B : Biomass (g/m2)
    B<-rep(NA,ldate)
    # LAI : Leaf Area Index (m2 leaf/m2 soil)
    LAI<-rep(NA,ldate)
    # CumInt LAI : Cumulative intercepted radiation
	CumInt<-rep(NA,ldate) 


    # Initialize state variables when sowing on day "sdate"   
    TT[sdate]<- 0
    B[sdate]<- 1
    LAI[sdate]<- 0.01
    CumInt[sdate] = 0.0

    # Simulation loop
    for (day in sdate:(ldate-1))
        {
        # Calculate rates of change of state variables (dTT, dB, dLAI)
        dTT <- max((weather$Tmin[day]+weather$Tmax[day])/2-Tbase, 0) 
		tday = (weather$Tmin[day]+weather$Tmax[day])/2
        if (TT[day]<=TTM) {dB <- maize.RUEtemp(tday,RUE_max,6.2,16.5,33,44)*(1-exp(-K*LAI[day]))*weather$I[day]}
        else {dB <- 0}
        if (TT[day]<=TTL) {dLAI <- alpha*dTT*LAI[day]*max(LAImax-LAI[day],0)}
        else {dLAI <-0 }
		
        
        # Uptade state variables 
        TT[day+1]<- TT[day] + dTT 
        B[day+1]<- B[day] + dB
        LAI[day+1]<- LAI[day] + dLAI   
		CumInt[day+1] = CumInt[day] + weather$I[day]*(1 - exp(- K * LAI[day]))
        }
        # End simulation loop
    return(data.frame(day=sdate:ldate,TT=TT[sdate:ldate],LAI=LAI[sdate:ldate],B=B[sdate:ldate],CumInt=CumInt[sdate:ldate]))    
    }

###############################################################################
#' @title The Maize model function with CumInt, RUE and ear growth
#' @param Tbase,RUE_max,K,alpha,LAImax,TTM,TTL : 7 parameters, weather table, sowing date and last date
#' @param  weather : weather data.frame for one single year
#' @param sdate : sowing date
#' @param ldate : last date
#' @return data.frame with daily TT, LAI,B
#' @export
maize_cir_rue_ear.model<-function(Tbase,RUE_max,K,alpha,LAImax,TTM,TTL,weather,sdate,ldate)
    {
    # Initialize variables
    # 3 states variables, as 3 vectors initialized to NA
    # TT : temperature sum (°C.d)
    TT<-rep(NA,ldate)
    # B : Biomass (g/m2)
    B<-rep(NA,ldate)
    # LAI : Leaf Area Index (m2 leaf/m2 soil)
    LAI<-rep(NA,ldate)
    # CumInt  Cumulative intercepted radiation
	CumInt<-rep(NA,ldate) 
	 # BE:  Biomass of ear
	BE<-rep(NA,ldate) 

    # Initialize state variables when sowing on day "sdate"
    TT[sdate]<- 0
    B[sdate]<- 1
    LAI[sdate]<- 0.01
    CumInt[sdate] = 0.0
	BE[sdate] = 0.0

    # Simulation loop
    for (day in sdate:(ldate-1))
        {
        # Calculate rates of change of state variables (dTT, dB, dLAI)
        dTT <- max((weather$Tmin[day]+weather$Tmax[day])/2-Tbase, 0) 
		tday = (weather$Tmin[day]+weather$Tmax[day])/2
        if (TT[day]<=TTM) {dB <- maize.RUEtemp(tday,RUE_max,6.2,16.5,33,44)*(1-exp(-K*LAI[day]))*weather$I[day]}
        else {dB <- 0}

        if (TT[day]<=TTL) {dLAI <- alpha*dTT*LAI[day]*max(LAImax-LAI[day],0)}
        else {dLAI <-0 }

	    if (TT[day]> TTL) {dBE <- dB}
		else dBE <-0			
		
        
        # Uptade state variables 
        TT[day+1]<- TT[day] + dTT 
        B[day+1]<- B[day] + dB
        LAI[day+1]<- LAI[day] + dLAI   
		CumInt[day+1] = CumInt[day] + weather$I[day]*(1 - exp(- K * LAI[day]))
		BE[day+1] <- BE[day] + dBE
        }
        # End simulation loop
    return(data.frame(day=sdate:ldate,TT=TT[sdate:ldate],LAI=LAI[sdate:ldate],B=B[sdate:ldate],CumInt=CumInt[sdate:ldate],BE=BE[sdate:ldate]))    
    }

###############################################################################
#' @title Reading Weather data function for MAIZE model (West of France Weather)
#' @param working.year : year for the subset of weather data (default=NA : all the year)
#' @param working.site : site for the subset of weather data (default=NA : all the site)
#' @param weather_all : weather data base (default=weather_FranceWest)
#' @return data.frame with daily weather data for one or several site(s) and for one or several year(s)
#' @export
# Reading Weather data function
maize.weather <- function(working.year=NA, working.site=NA,weather_all=weather_FranceWest)
    {
    # WEYR => year
    # WEDAY => day
    # SRAD => I: solar radiation (MJ)
    # TMAX => Tmax : maximum temperature (°C)
    # TMIN => Tmin : minimum temperature (°C)
    # select only useful header
    weather=weather_all[,c("idsite","GPSlatitude","GPSlongitude","WEYR","WEDAY","SRAD","TMAX","TMIN")]
    names(weather)[names(weather)=="WEDAY"]= "day"
    names(weather)[names(weather)=="WEYR"]= "year"
    names(weather)[names(weather)=="SRAD"]= "I"
    names(weather)[names(weather)=="TMAX"]= "Tmax"
    names(weather)[names(weather)=="TMIN"]= "Tmin"
    # if argument working.year/working.site is specified, work on one particular year/site
    if (!is.na(working.year)&!is.na(working.site)) {weather=weather[(weather$year==working.year)&(weather$idsite==working.site),] }
    else{
      if (!is.na(working.year)) {weather=weather[(weather$year==working.year),]}
      if (!is.na(working.site)) {weather=weather[(weather$idsite==working.site),]}}
    return (weather)
    }
################################################################################
# End of file