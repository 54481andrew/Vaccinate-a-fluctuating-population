## This script outputs a figure that shows how the number of susceptible,
## infected, and recovered hosts varies over the course of a year when the
## population has reached a stable limit cycle. This script prints the
## size of each class at time t%%T = 0, which is then fed into the Gillespie
## algorithm as an initial condition. 

rm(list = ls(all = TRUE))

############################################
############# USER INPUT ###################
############################################

parmat <- expand.grid(d = c(1/(1*365)),
                      gamv = 1/14,
                      Nv = 0,
                      tb = 120, T = 365, NPeak = 2696,
                      gamp = c(1/(30)), pmu = 0,
		      Rp = c(2),
                      tv = 0, type = c('f'))

## Use the specified peak population size to
## determine the constant birth rate during
## the breeding season. This relationship
## NPeak and b is derived in the supplemental
## Mathematical file.
parmat$b     <- with(parmat,
                     NPeak*d/(exp(d*(T-tb))*(exp(d*tb)-1)/(exp(d*T)-1)))
## The transmission coefficient depends on the mode of transmission that is
## assumed. 'f' is frequency dependent, 'd' is density dependent. 
parmat$Bp[parmat$type=='f']    <- with(parmat[parmat$type=='f',], Rp*(d+gamp))
parmat$Bp[parmat$type=='d']    <- with(parmat[parmat$type=='d',], Rp*d*(d+gamp)*T/(b*tb))

############################################
############## END USER INPUT ##############
############################################

NPathVals <- nrow(parmat)

## Specify the amount of time that the system is simulated - 
## should be long enough to ensure stable limit cycle is
## reached. 
TBurnIn <- 10*365

## Load in a package that allows parallel computation, and
## the solving of ode's with lsoda. Also, load in functions
## that return the differentials of the system of ode's that
## our paper studies. 
require(parallel)
require(deSolve)
source('Tools/Functions.r', local = TRUE)


#########################################
#### DEFINE THE SIMULATION FUNCTION #####
#########################################

## This function performs the simulation.
## mcl.fun takes an argument "i" that tells the
## simulation to use parameters from the i^th
## row of the parameter matrix.

starttime <- Sys.time()
mclBurn <- function(i){
    
    ## Define a vector of times at which vaccination occurs
    tv <- parmat$tv[i]
    T  <- parmat$T[i]
    burntimes <- seq(0, TBurnIn, by = 1)
    Sinit  <- 1000
    Svinit <- 0
    Ipinit <- 100
    Vinit  <- 0
    Pinit  <- 0
    yinit  <- c(Sinit, Svinit, Ipinit, Vinit, Pinit)

    rhs <- if(parmat$type[i]=='f'){rhs.freq.full}else{rhs.dens.full}
    out       <- lsoda(y = yinit, times = burntimes,
                       hmax = 0.1, func = rhs,
                       parms = parmat[i,], maxsteps = 1000000)
    return(out)
}

BurnSolSet <- mclapply(X = 1:NPathVals, FUN = mclBurn, mc.cores = detectCores() - 2)

## Graph
tc <- 50
wi.graph <- 1
for(i in wi.graph){
    dval <- parmat$d[i]
    Rpval <- parmat$Rp[i]
    out <- data.frame(BurnSolSet[[i]])
    names(out) <- c('time', 'S', 'Iv', 'Ip', 'V', 'P')
    NPeak = parmat$NPeak[i]
    out$N = with(out, S + Iv + Ip + V + P)
    mod.out <- out[with(out, time >= TBurnIn-365),]
    mod.out <- mod.out[seq(0,nrow(mod.out),length.out = 25),]
    mod.out$N = with(mod.out, S + Iv + Ip + V + P)
    mod.out$pIp = with(mod.out, Ip/N)
    mod.out$pTotInf = with(mod.out, (Ip + P)/N)
    mod.out$pS = with(mod.out, S/N)
    mod.out <- mod.out[-nrow(mod.out),]
    mod.out$time <- mod.out$time%%365
    ## Population Size
    setEPS()
    fn <- 'Figure_Rat.eps'
    postscript(file = fn, width = 5,
        height = 4)
    par(mai = c(0.75, 0.75, 0.25, 0.25))
    plot(NA, xlim = c(0,365), ylim = c(0,1.2), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    tb <- unique(parmat$tb)
    polygon(x = c(0,tb,tb,0), y =c(-100,-100,1.2*1.03,1.2*1.03), col = 'gray90',
            border = NA)
    matlines(mod.out$time, mod.out$Ip/NPeak, type = 'b', pch = 2,
             lty = 1, lwd = 2, col = 'red', cex = 1)
    matlines(mod.out$time, (mod.out$Ip + mod.out$P)/NPeak, type = 'b', pch = 2,
        lty = 1, lwd = 2, col = 'gray', cex = 1)
    matlines(mod.out$time, mod.out$S/NPeak, type = 'b', pch = 2,
             lty = 1, lwd = 2, col = 'blue', cex = 1)
    wi.lastyear <-  max(out$time) - 365 < out$time  & out$time < max(out$time)
    matlines(mod.out$time, mod.out$N/NPeak, lwd = 2)
    abline(v = seq(0,365, by = 30.4), lwd = 0.25, lty = 3)
    abline(h = seq(0,1.2, by = 0.1), lwd = 0.25, lty = 3)
    legend(x = 'topleft',legend = c('Total Pop Size', 'Susceptible', 'Lassa Infected', 'Recovered'),
           col = c('black', 'blue', 'red', 'gray'), pch = c(NA,2,2,2,1,1), lwd = 1.5, cex = 0.8, bg = 'white', ncol = 2)
    mtext(side = 1, text = '', line = 2.5)
mtext(side = 2, text = 'Fraction of Peak Population', line = 2.5)
    axis(side = 1, labels = FALSE, at = seq(0,365, by = 365/12) )
    axis(side = 2, labels = TRUE, at = seq(0,1, by = 0.2))
    text(x = -10 + seq(0,365, length.out = 13), -0.1, labels = month.abb[c(6:12,1:6)],
         srt = 65, pos = 1, xpd = TRUE)
    dev.off()
}## Close figure loop

## Calculate and print the size of each class at t%%T = 0. 
wi <- (out$time == TBurnIn)
state <- round(out[wi,])
state$N <- with(state, S + Iv + Ip + V + P)
print(state)


