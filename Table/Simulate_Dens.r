## This script calculates the vaccination window for the species listed in
## Table 3 of the main text.

rm(list = ls(all = TRUE))

## Load package that allows for parallel computation, numerical
## simulation of ODE's with lsoda, and additional functions that
## compute the differentials of the ODE system.
require(parallel)
require(deSolve)
source('Tools/Functions.r', local = TRUE)

datname = "Animals_Density" ## Name of the data file that is written

## Load in a parameter matrix that contains demographic
## parameters on various wildlife.
source("make_parmat.r", local = TRUE)

## Add in row names
NPars = nrow(parmat)
parmat$NPar = 1:NPars

## Use the specified peak population size to
## determine the constant birth rate during
## the breeding season. This relationship
## NPeak and b is derived in the supplemental
## Mathematical file.
parmat$b     <- with(parmat,
                     NPeak*d/(exp(d*(T-tb))*(exp(d*tb)-1)/(exp(d*T)-1)))


nyears   <- 100 ## Number of years to simulate each parameter set
maxtimes <- 365*nyears
times    <- seq(0,maxtimes, by = 1)
witimes <- times >= (nyears-1)*365 ## Times used to assess the vaccination program
yeartimes <- c(times[witimes]%%365)
yeartimes[length(yeartimes)] <- 365
tvvals <- matrix(seq(0.1,364.9,length.out = 50), ncol = 1) ## Used to approximate V vs tv fxn (below)

## Write the parameter matrix
filename = paste("Data/", datname, sep = '')
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

    y0 = c(1000, 0, 0)
    names(y0) = c('S','Sv','V')
    TVal  <- parmat$T[i]

    cost.fxn <- function(tv,parms){
        NVacc <- parms$Nv
        vacctimes <- seq(tv, maxtimes,by = TVal)

        out  <- data.frame(lsoda(y = y0, times = times,
                                 func = rhs.jv, parms = parms,
                                 events=list(func = vaccinate.jv, time=vacctimes),
                                 maxsteps = 100000))

        ## Return average serorpevalence
	vals    <- with(out[witimes,], V)
        vals[1] <- vals[length(vals)]
	vapprox = splinefun(x = yeartimes, vals, method = "periodic")
        npop <- with(out[witimes,], (S + Sv + V))
        npopapprox = splinefun(x = yeartimes, npop, method = "periodic")
        VAvg <- 1/365*(integrate(vapprox, lower = 0, upper = 365, stop.on.error = FALSE)$value)
        NAvg <- 1/365*(integrate(npopapprox, lower = 0, upper = 365, stop.on.error = FALSE)$value)
        return(VAvg/NAvg)
    }## End cost.fxn


    ## In order to calculate the vaccination window, we must know how well
    ## suboptimal vaccination programs perform. Here, we evaluate cost.fun
    ## across a range of tv values. Later, this data will be used to guide
    ## the optimize function to find the best vaccination strategy.
    vnullvals <- apply(X = tvvals, MARGIN = 1, FUN = cost.fxn, parms = parmat[i,])
    tvstar <- tvvals[which.max(vnullvals)]

    ## Use vnullvals to find the interval that contains the best strategy (maximum VAvg)
    wi.interval <- c(max(1, which.max(vnullvals) - 1), min(which.max(vnullvals) + 1, length(tvvals)))
    interval <- tvvals[wi.interval]

    ## Find optimal time of vaccination and corresponding average seroprevalence
    uroot <- optimize(f=cost.fxn, interval = interval, parmat[i,], maximum=TRUE)
    tvstar <- uroot$maximum ## Best tv value
    vavgstar <- uroot$objective ## Maximum average seroprevalence

    ## Find and store the 90% and 95% vaccination window
    wi90 <- which(vnullvals/vavgstar > 0.9)
    parmat$downtv90[i] <- min(tvvals[wi90])
    parmat$uptv90[i] <- max(tvvals[wi90])

    wi95 <- which(vnullvals/vavgstar > 0.95)
    parmat$downtv95[i] <- min(tvvals[wi95])
    parmat$uptv95[i] <- max(tvvals[wi95])

    ## Store the optimal strategy (tv, vavg) in parmat
    parmat$tvstar[i] <- tvstar
    parmat$vavgstar[i] <- vavgstar

    ## Store the typical average seroprevalence across all vaccinations in parmat
    vnullfun = splinefun(x = c(tvvals,365.1), c(vnullvals,vnullvals[1]))
    parmat$vnull[i]  <- 1/365*integrate(vnullfun, lower = 0, upper = 365)$value

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

refparmat <- read.table(file = filename, header = T)
refparmat <- refparmat[order(refparmat$NPar),]
write.table(refparmat, file = filename, append = F, row.names = F)
parmat <- read.table(file = filename, header = T)

parmat$Window90 <- with(parmat, uptv90 - downtv90)
parmat$Window95 <- with(parmat, uptv95 - downtv95)

write.table(parmat, file = filename, append = F, row.names = F)

tab <- subset(parmat, select = c('Animal', 'Window95'))
tab$Month95 <- round(tab$Window95/30.4,1)
tab[,c(1,3)]
