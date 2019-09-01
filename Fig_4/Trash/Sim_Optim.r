## This script calculates the fractional reduction in mean abundance
## of pathogen infected hosts.


############################################
############# USER INPUT ###################
############################################


rm(list = ls(all = TRUE))

SimName <- 'Fig_4'

parmat <- expand.grid(d = c(1/(1*365)),
                      gamv = c(1/14, 1/30),
                      Nv = c(250, 500),
                      tb = 120, T = 365, NPeak = 1000,
                      gamp = c(1/30.4, 1/365), pmu = c(0,0.9),
		      Rp = c(2,5),
                      tv = seq(0, 365, length.out = 25),
                      type = c('f', 'd'))

## Use the specified peak population size to
## determine the constant birth rate during
## the breeding season. This relationship
## NPeak and b is derived in the supplemental
## Mathematical file.
parmat$b     <- with(parmat,
                NPeak*d/(exp(d*(T-tb))*(exp(d*tb)-1)/(exp(d*T)-1)))
parmat$NAvg  <- with( parmat, b*tb/(d*T) ) ## Average population size in absence of pathogen

## Calculation of Bp from Rp depends on mode of transmission
parmat$Bp    <- with(parmat, ifelse(type=='f', Rp*(d+gamp), Rp*(d+gamp)/NAvg ))

############################################
############## END USER INPUT ##############
############################################

## Load packages for lsoda and parallel
## simulations. Functions.r for functions that
## return differentials of the ODE.
require(deSolve)
require(parallel)
source('Tools/Functions.r', local = TRUE)

## Give each row of parmat a name
NPars = nrow(parmat)
parmat$NPar = 1:NPars

## Add columns that store data from the simulation
parmat$Ip0 <- NA
parmat$Ipv <- NA


nyears   <- 100
TBurnIn  <- 365*95
maxtimes <- 365*nyears
times    <- seq(0,maxtimes, by = 1)


witimes <- times >= TBurnIn ## begin assessing the campaign on year 95
yeartimes <- c(times[witimes]%%365)
yeartimes[length(yeartimes)] <- 365


###############################################
########## SET UP DATA STORAGE ################
###############################################

## Create a text file so that for each parameter
## set, the parameter values, simulation errors,
## and simulation output are saved. The resulting
## data file is stored within the "Data" folder.

filename = paste("Data/", datname, sep = '')
## Write the parameter matrix
write.table(parmat[1,][-1,], file = filename, append = F)


#########################################
#### DEFINE THE SIMULATION FUNCTION #####
#########################################

## This function simulates the ODE system, and uses
## the "optimize" function to find the timing of
## vaccination that best elevates the population's
## seroprevalence.
## mcl.fun takes an argument "i" that tells the
## simulation to use parameters from the i^th
## row of the parameter matrix.

mcl.fun = function(i){


    tvval <- parmat$tv[i]


    NPeak <- parmat$NPeak[i]
    ## Initial condition. Begin simulation with susceptible population. Add in 5 infected hosts.
    y0 = with(parmat[i,], c(NPeak*exp(-d*(365-tb)), 0, 5, 0, 0))
    names(y0) = c('S','Iv','Ip','V','P')

    NVacc          <- parmat$Nv[i]
    parmat$Nv[i] <- 0 #temporarily set to 0

    source("Tools/Functions.r", local = T)
    rhs <- if( parmat$type[i] == 'f'){ rhs.freq.full }else{ rhs.dens.full }

    out1 <- data.frame(lsoda(y = y0, times = seq(0,TBurnIn,by=1),
                             func = rhs, parms = parmat[i,]))
    y0 = out1[nrow(out1), -1]
    y0 = as.numeric(y0)
    parmat$Nv[i] <- NVacc
    names(y0) = c('S','Iv','Ip','V','P')
    vacctimes <- seq(TBurnIn + tvval, maxtimes + tvval + 50*365, by = 365)
    out2 <- data.frame(lsoda(y = y0,
                             times = seq(TBurnIn, maxtimes + tvval + 50*365,
                                                 by = 1),
                             func = rhs, parms = parmat[i,],
                             events=list(func = vaccinate.full, time=vacctimes),
                             maxsteps = 100000))

    ## Extract information to be graphed

    ## Extract mean pathogen abundance in the last 5 years of TBurnIn period,
    ## and 5 years after.
    Ip.novacc.mean <- mean(out1$Ip[out1$time >= TBurnIn - 5*365])

    ## Method 1: fixed 5 year interval
    Ip.withvacc.mean <- mean( out2$Ip[out2$time >= TBurnIn &
                                      out2$time <= TBurnIn + 5*365] )
    ##Method 2: rolling 5 year interval
    Ip.withvacc.mean.alt <- mean( out2$Ip[out2$time >= (TBurnIn + tvval) &
                                          out2$time <= TBurnIn + 5*365 + tvval] )
    ## Method 3: longterm avg - 50 years
    Ip.withvacc.mean.long <- mean( out2$Ip[out2$time >= TBurnIn] )

    ## store into parmat
    parmat[i,'Ip0'] <- Ip.novacc.mean
    parmat[i,'Ipv'] <- Ip.withvacc.mean
    parmat[i,'Ipvalt'] <- Ip.withvacc.mean.alt
    parmat[i,'Ipvlong'] <- Ip.withvacc.mean.long

    ## write data
    write.table(parmat[i,], file = filename, append = T, col.names = F, row.names = F)
    print(paste("NPar:", i, ' / ', NPars, ' ; Time: ', Sys.time() - starttime),
          quote = FALSE)

    return(i)
}

starttime <- Sys.time()
mcl.out <- mclapply(X = 1:NPars, FUN = mcl.fun, mc.cores = detectCores()-1)

refparmat <- read.table(file = filename, header = T)
refparmat <- refparmat[order(refparmat$NPar),]
write.table(refparmat, file = filename, append = F, row.names = F)
parmat <- read.table(file = filename, header = T)

