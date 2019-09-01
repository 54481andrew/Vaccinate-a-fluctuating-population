## Define a function that generates partially transparent colors
## color = color name
## percent = % transparency
## name = an optional name for the color
t_col <- function(color, percent = 50, name = NULL) {
    ## Get RGB values for named color
    rgb.val <- col2rgb(color)
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100-percent)*255/100,
                 names = name)
    ## Save the color
    invisible(t.col)    
}

## Define a function that returns the differentials of the
## ODE system, assuming that transmission is density dependent
rhs.dens.full = function(t,y,parms){
    S <- y[1]
    Sv <- y[2]
    Ip <- y[3]
    V <- y[4]
    P <- y[5]
    with(parms,{
        dS = b*(t%%T < tb) - Bp*S*Ip - d*S
        dSv = -(gamv + d)*Sv - Bp*Sv*Ip
        dIp = Bp*S*Ip + Bp*Sv*Ip - (gamp + d)*Ip
        dV = gamv*Sv - d*V
        dP = gamp*(1 - pmu)*Ip - d*P
        
        return(list(c(dS,dSv,dIp,dV,dP)))
    })
}

## Define a function that returns the differentials of the
## ODE system, assuming that transmission is frequency dependent.
rhs.freq.full = function(t,y,parms){
    S <- y[1]
    Sv <- y[2]
    Ip <- y[3]
    V <- y[4]
    P <- y[5]
    NPop <- sum(y)
    with(parms,{
        dS = b*(t%%T < tb) - Bp*S*Ip/NPop - d*S
        dSv = -(gamv + d)*Sv - Bp*Sv*Ip/NPop
        dIp = (Bp*S*Ip + Bp*Sv*Ip)/NPop - (gamp + d)*Ip
        dV = gamv*Sv - d*V
        dP = gamp*(1-pmu)*Ip - d*P        
        return(list(c(dS,dSv,dIp,dV,dP)))
    })
}

## Define the function that is called during
## an event in lsoda. This function modifies
## the state variables in accordance with a
## single pulse vaccination.
vaccinate.full <- function(t,y,parms){
    S <- y[1]
    Sv <- y[2]
    Ip <- y[3]
    V <- y[4]
    P <- y[5]
    NPop <- sum(y)
    with(parms,{
	nvacc <- min(S, Nv*S/NPop)
	S   <- S - nvacc
	Sv  <- Sv + nvacc
	return(c(S, Sv, Ip, V, P))
	})
}


## Define a function that returns the differentials
## of the ODE system in the absence of the pathogen
## (Ip = 0, P = 0). 
rhs.jv = function(t,y,parms){
    S <- y[1]
    Sv <- y[2]
    V <- y[3]
    with(parms,{
        dS = b*(t%%T < tb) - d*S
        dSv = -(gamv + d)*Sv
        dV = gamv*Sv - d*V        
        return(list(c(dS,dSv,dV)))
    })
}

## Define the function that is called during lsoda
## when the pathogen is absent. Similar to vaccinate.full,
## but Ip and P are assumed 0. 
vaccinate.jv = function(t,y,parms){
    S <- y[1]
    Sv <- y[2]
    V <- y[3]
    NPop <- sum(y)
    with(parms,{
        nvacc <- round(min(S, Nv*S/NPop))
	S   <- S - nvacc
	Sv  <- Sv + nvacc
	return(c(S, Sv, V))
   })
}
