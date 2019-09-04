## This script calculates the fractional reduction in mean abundance
## of pathogen infected hosts.


############################################
############# USER INPUT ###################
############################################


rm(list = ls(all = TRUE))
datname <- 'Fig_4'
parmat <- expand.grid(d = c(1/(1*365)),
                      gamv = c(1/14, 1/30),
                      Nv = c(250),
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

SimTime.path  <- 365*50 ## Years taken for pathogen to reach stable cycle
SimTime.vacc  <- 365*50 ## Years taken for vaccination to reach stable cycle

times.path.sim <- seq(0,SimTime.path - 1, by = 1) ## times argument for lsoda, pathogen only simulation
times.vacc.sim <- seq(SimTime.path, SimTime.path + SimTime.vacc, by = 1) ## times argument for vaccination simulation


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

    ## Vaccination time used for this simulation
    tvval <- parmat$tv[i]

    ## Initial condition. Begin simulation with susceptible population. Add in 5 infected hosts.
    y0 = with(parmat[i,], c(NPeak*exp(-d*(365-tb)), 0, 5, 0, 0))
    names(y0) = c('S','Iv','Ip','V','P')

    ## Temporarily set vaccines to 0 for pathogen-only simulation
    NVacc <- parmat$Nv[i]
    parmat$Nv[i] <- 0

    ## Load in the correct functional form for frequency or density dependent transmission
    source("Tools/Functions.r", local = T)
    rhs <- if( parmat$type[i] == 'f'){ rhs.freq.full }else{ rhs.dens.full }

    ## Simulate the system with S, Ip, and P classes only
    out1 <- data.frame(lsoda(y = y0, times = times.path.sim,
                             func = rhs, parms = parmat[i,]))
    ## Reset the initial condition to the state vector at the end of the pathogen-only simulation
    y0 = out1[nrow(out1), -1]
    y0 = as.numeric(y0)
    names(y0) = c('S','Iv','Ip','V','P')

    ## Add vaccines back in, and specify the vaccination times
    parmat$Nv[i] <- NVacc
    vacctimes <- seq(SimTime.path + tvval, max(times.vacc.sim), by = 365)

    ## Simulate forward with vaccination
    out2 <- data.frame(lsoda(y = y0,
                             times = times.vacc.sim,
                             func = rhs, parms = parmat[i,],
                             events=list(func = vaccinate.full, time=vacctimes),
                             maxsteps = 100000))

    ## Extract information to be graphed

    ## Extract mean pathogen abundance in the last 25 years of SimTime period,
    ## and 25 years following 25 years of vaccination.
    Ip.novacc.mean <- mean(out1$Ip[out1$time >= SimTime.path - 25*365])
    Ip.withvacc.mean <- mean( out2$Ip[out2$time >= SimTime.vacc + 25*365] )

    ## Store into parmat
    parmat[i,'Ip0'] <- Ip.novacc.mean
    parmat[i,'Ipv'] <- Ip.withvacc.mean

    ## Write the data by adding a the data from row i of parmat
    ## to the data file [filename]
    write.table(parmat[i,], file = filename, append = T, col.names = F, row.names = F)
    print(paste("NPar:", i, ' / ', NPars, ' ; Time: ', Sys.time() - starttime),
          quote = FALSE)
    return(i)
}

## Run the simulations. If mclapply throws an error, lapply can be used.
starttime <- Sys.time()
ncores <- detectCores()-2 ## Number of cores to use
mcl.out <- mclapply(X = 1:NPars, FUN = mcl.fun, mc.cores = ncores)
tottime <- Sys.time() - starttime
print(paste('Total Time: ', tottime, ' ; ' , 'CPU time: ', ncores*tottime, sep = ''), quote = FALSE)


## Read in the simulated data, reorder the rows, and rewrite to file
refparmat <- read.table(file = filename, header = T)
refparmat <- refparmat[order(refparmat$NPar),]
write.table(refparmat, file = filename, append = F, row.names = F)
parmat <- read.table(file = filename, header = T)

