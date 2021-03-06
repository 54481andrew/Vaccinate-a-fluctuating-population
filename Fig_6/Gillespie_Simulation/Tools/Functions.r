### Function to generate transparent colors
t_col <- function(color, percent = 50, name = NULL) {
                                        #	  color = color name
                                        #	percent = % transparency
                                        #	   name = an optional name for the color
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
        dP = gamp*Ip - d*P

        return(list(c(dS,dSv,dIp,dV,dP)))
    })
}

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
        dP = gamp*Ip - d*P
        return(list(c(dS,dSv,dIp,dV,dP)))
    })
}

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


###Analogous functions with just vaccination (jv)

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
