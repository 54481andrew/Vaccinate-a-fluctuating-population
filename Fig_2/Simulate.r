## This code simulates a pulse vaccination
## program.

############################################
############# USER INPUT ###################
############################################


rm(list = ls(all=TRUE))
datname  <- "Fig_2" ## Name of the data file produced
parmat <- expand.grid(d = c(1/(365*1),
                            1/(365*2.5), 1/(365*5), 1/(365*10)),
                      gamv = 1/14,
                      Nv = 500,
                      tv = seq(0,365,length.out = 100),
                      tb = c(30, 90, 180, 270), T = 365, NPeak = 1000)

## Use the specified peak population size to
## determine the constant birth rate during
## the breeding season. This relationship
## NPeak and b is derived in the supplemental
## Mathematical file.
parmat$b     <- with(parmat, NPeak*d/(exp(d*(T-tb))*(exp(d*tb)-1)/(exp(d*T)-1)))



############################################
############## END USER INPUT ##############
############################################

## Load packages for lsoda and parallel
## simulations.
require(deSolve)
require(parallel)

###Add columns that store the row number,
###and statistics on the min, average, and
###max of S and V
npars <- nrow(parmat)
parmat$npar <- 1:npars
parmat$MaxSusfrac <- NA
parmat$MaxVfrac <- NA
parmat$AvgSusfrac <- NA
parmat$AvgVfrac <- NA
parmat$MinSusfrac <- NA
parmat$MinVfrac <- NA
parmat$VAbund <- NA
parmat$NPopAbund <- NA
parmat$VdivNPop <- NA

maxtimes <- 100*365 #Years
ftimeseq <- seq(maxtimes-365,maxtimes,by = 1)
timeseq  <- c(0,maxtimes-2*365, ftimeseq )
tseq     <- ftimeseq - (maxtimes-365)

###############################################
########## SET UP DATA STORAGE ################
###############################################

## Create a text file so that for each parameter
## set, the parameter values, simulation errors,
## and simulation output are saved. The resulting
## data file is stored within the "Data" folder.

filedatname <- paste("Data/",datname,sep="")
## Write the parameter matrix
write.table(parmat[1,][-1,], file = filedatname, append = FALSE,
            row.names = FALSE, col.names = TRUE)


#########################################
#### DEFINE THE SIMULATION FUNCTION #####
#########################################

## This function simulates the ODE system.
## mcl.fun takes an argument "i" that tells the
## simulation to use parameters from the i^th
## row of the parameter matrix.

mcl.fun <- function(i){

    ## Set initial condition
    Sinit  <- 1000
    Svinit <- 0
    Vinit  <- 0
    yinit  <- c(Sinit, Svinit, Vinit)

    ## Load functions that return differentials
    source('Tools/Functions.r', local = TRUE)

    ## Define a vector of times at which vaccination occurs
    tv <- parmat$tv[i]
    T  <- parmat$T[i]
    vacctimes <- seq(tv, maxtimes, by = T)

    ## lsoda simulates our ODE (defined by rhs.jv). Note
    ## that in addition to the argument "times", lsoda also
    ## has an argument "events". This tells lsoda to stop
    ## the simulation and apply the vaccinate.jv function
    ## to the state variables.
    out       <- lsoda(y = yinit, times = timeseq,
                       hmax = 0.1, func = rhs.jv,
                       parms = parmat[i,],
                       events=list(func = vaccinate.jv, time = vacctimes))

    ## Get the solution at the stable cycle
    ind <- out[,1] %in% ftimeseq
    y   <- out[ind,-1]
    Svals  <- y[,1]
    Svvals <- y[,2]
    Vvals  <- y[,3]

    ## Store the solution as a dataframe, add columns to the
    ## data frame that give the fraction of seronative and seropositve
    ## individuals.
    sol <- data.frame(day = ftimeseq, S = Svals, Sv = Svvals,
                      V = Vvals, NPop = Svals + Svvals + Vvals)
    sol$FracVacc <- with(sol, V/NPop)
    sol$FracSusc <- with(sol, 1 - FracVacc)

    ## Get statistics from simulation using a periodic spline interpolants.
    ## The interpolant is fit to the final year of simulation time.
    ind  <- which(sol$day >= maxtimes-365)
    SeroNegfun <- splinefun(sol[ind,c('day','FracSusc')], method = "periodic")
    Vfun <- splinefun(sol[ind,c('day','V')], method = "periodic")
    NPopfun <- splinefun(sol[ind,c('day','NPop')], method = "periodic")

    parmat$MinSusfrac[i] <- min(SeroNegfun(tseq))
    parmat$MaxVfrac[i] <- 1 - min(SeroNegfun(tseq))
    parmat$MaxSusfrac[i] <- max(SeroNegfun(tseq))
    parmat$MinVfrac[i] <- 1 - max(SeroNegfun(tseq))
    parmat$AvgSusfrac[i] <- mean(SeroNegfun(tseq))
    parmat$AvgVfrac[i] <- 1 - mean(SeroNegfun(tseq))
    parmat$VAbund[i] <- mean(Vfun(tseq))
    parmat$NPopAbund[i] <- mean(NPopfun(tseq))
    parmat$VdivNPop[i] <- mean(Vfun(tseq))/mean(NPopfun(tseq))

    ## Write the data by adding a
    ## the data from row i of parmat
    ## to the data file [filedatname]
    write.table(parmat[i,], file = filedatname,
                append = TRUE, row.names = FALSE, col.names = FALSE)

    ## Print and return index to help with debugging
    print(paste("npar: ", i, " / ", npars, ';  Time: ',Sys.time() - starttime, sep=""))
    return(i)

} ##Close the function mcl.fun

## Simulate the system. If mclapply throws an error, try
## using the lapply code that is commented out.
starttime <- Sys.time()
ncores <- detectCores()-2 ## Number of cores to use
mcl.out <- mclapply(X = 1:NPars, FUN = mcl.fun, mc.cores = ncores)
## mcl.out <- lapply(X=1, FUN = mcl.fun)
tottime <- Sys.time() - starttime
print(paste('Total Time: ', tottime, ' ; ' , 'CPU time: ', ncores*tottime, sep = ''), quote = FALSE)

## Read in the simulated data, reorder the rows, and rewrite to file
parmat_filled <- read.table(file = filedatname, header = TRUE)
parmat_filled_ordered <- parmat_filled[order(parmat_filled$npar),]
write.table(parmat_filled_ordered, file = filedatname,
            append = FALSE, row.names = FALSE, col.names = TRUE)

