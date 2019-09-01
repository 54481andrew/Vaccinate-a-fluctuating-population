############################################
###This code simulates a pulse vaccination
###program. The first pulse occurs at time VaccT0,
###and subsequent pulses occur every VaccPer years.
###The population birth rate is b between times
###TB1 and TB2, and 0 elsewhere.
#*******************************************

### Each vaccination vaccinates a fixed
### number of individuals. The number of
### vaccinations is NVacc. In parmat, this
### is also described as rho, the ratio of
### NVacc to the peak population size (NPeak).

############################################
############# USER INPUT ###################
############################################

rm(list = ls(all=TRUE))
## Simulation time: 74 minutes cpu time

datname  <- "Fig_3" ###Name of the data file
parmat <- expand.grid(d = c(1/(1*365)),
                      gamv = 1/14,
                      rho = c(0.25, 0.5, 0.75),
                      tv = seq(0,365,length.out = 25),
                      tb = c(3*30.4), T = 365, NPeak = 30000)
parmat$b     <- with(parmat, NPeak*d/(exp(d*(T-tb))*(exp(d*tb)-1)/(exp(d*T)-1)))
parmat$Nv    <- with(parmat, NPeak*rho)

############################################
############## END USER INPUT ##############
############################################
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

maxtimes <- 40*365 #Years
ftimeseq <- seq(maxtimes-365,maxtimes,by = 1)
timeseq  <- c(0,maxtimes-2*365, ftimeseq )
tseq     <- ftimeseq - (maxtimes-365)


###############################################
########## SET UP DATA STORAGE ################
###############################################
###Create a text file so that, for each parameter
###set, the parameter values, simulation errors,
###and simulation output are saved. The resulting
###data file is stored within the "Data" folder.
###############################################
filedatname <- paste("Data/",datname,sep="")
wflag <- TRUE
if(wflag){
### Write the parameter matrix
    write.table(parmat[1,][-1,], file = filedatname, append = FALSE,
                row.names = FALSE, col.names = TRUE)
### dir.create(filematname)
}

#########################################
#### DEFINE THE SIMULATION FUNCTION #####
#########################################
###This function does all the simulation.
###mcl.fun takes an argument "i" that tells the
###simulation to use parameters from the i^th
###row of the parameter matrix.
mcl.fun <- function(i){

###Set initial condition as first vaccination
    Sinit  <- 100
    Svinit  <- 0
    Vinit  <- 0
    yinit  <- c(Sinit, Svinit, Vinit)

###Define functions that return differentials
    source('Tools/Functions.r', local = TRUE)

###Define a vector of times at which vaccination occurs
    tv <- parmat$tv[i]
    T  <- parmat$T[i]
    vacctimes <- seq(tv, maxtimes, by = T)

###lsoda simulates our ODE (defined by rhs.ode). Note
###that in addition to the argument "times", lsoda also
###has an argument "events". This tells lsoda to stop
###simulation at each element of vacctimes, vaccinate,
###and continue simulation.
    out       <- lsoda(y = yinit, times = timeseq,
                           hmax = 0.1, func = rhs.jv,
                           parms = parmat[i,],
                           events=list(func = vaccinate.jv, time = vacctimes))
###Get the solution at the stable cycle
    ind <- out[,1] %in% ftimeseq
    y   <- out[ind,-1]
    Svals <- y[,1]
    Svvals <- y[,2]
    Vvals <- y[,3]

###Store the solution as a dataframe, add columns to the
###data frame that give the fraction of seronative and seropositve
###individuals.
    sol <- data.frame(day = ftimeseq, S = Svals, Sv = Svvals,
                      V = Vvals, Npop = Svals + Svvals + Vvals)
    sol$NPop     <- with(sol, S + Sv + V)
    sol$FracVacc <- with(sol, V/NPop)
    sol$FracSusc <- with(sol, 1 - FracVacc)

	###Get statistics from simulation
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

      ###Write the data by adding a
      ###the data from row i of parmat
      ###to the data file [filedatname]
      if(wflag){
          write.table(parmat[i,], file = filedatname,
                    append = TRUE, row.names = FALSE, col.names = FALSE)
      }


###Return something to help with debugging
    print(paste("npar: ", i, " / ", npars, ';  Time: ',Sys.time() - starttime, sep=""))
    return(i)
}###Close the function mcl.fun


###Here's where we tell R to simulate the system. If mclapply
###doesn't work, try using the lapply code that is commented out.
starttime <- Sys.time()
sollist <- mclapply(X=1:npars, FUN=mcl.fun, mc.cores=detectCores() - 2)
###sollist <- lapply(X=1, FUN = mcl.fun)
print(paste('Total Time: ', Sys.time() - starttime))

parmat_filled <- read.table(file = filedatname, header = TRUE)
parmat_filled_ordered <- parmat_filled[order(parmat_filled$npar),]
write.table(parmat_filled_ordered, file = filedatname,
            append = FALSE, row.names = FALSE, col.names = TRUE)



