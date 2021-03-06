## This code simulates a pulse vaccination
## program. 

############################################
############# USER INPUT ###################
############################################


rm(list = ls(all=TRUE))
datname  <- "Fig_1" ###Name of the data file produced
parmat <- expand.grid(d = c(1/(2.5*365)),
                      gamv = 1/14,
                      Nv = 500, 
                      tv = c(1*30.4, 5*30.4, 9*30.4),
                      tb = c(2*30.4, 12*30.4), T = 365, NPeak = 1000)
## Use the specified peak population size to
## determine the constant birth rate during
## the breeding season. This relationship
## NPeak and b is derived in the supplemental
## Mathematical file. 
parmat$b     <- with(parmat,
	     	(d*exp(-d*(T-tb))*(-1+exp(d*T))*NPeak)/(-1+exp(d*tb)))

############################################
############## END USER INPUT ##############
############################################

## Load packages for lsoda and parallel
## simulations. 
require(deSolve)
require(parallel)

## Add columns that store the row number,
## and statistics on the min, average, and
## max of S and V
npars <- nrow(parmat)
parmat$npar <- 1:npars
parmat$MaxSusfrac <- NA
parmat$MaxVfrac <- NA
parmat$AvgSusfrac <- NA
parmat$AvgVfrac <- NA
parmat$MinSusfrac <- NA
parmat$MinVfrac <- NA

maxtimes <- 100*365 # Number of years to simulate each parameter set
ftimeseq <- seq(maxtimes-365,maxtimes,by = 1) ## Times that are plotted in the figure
timeseq  <- c(0,maxtimes-2*365, ftimeseq ) ## "times" argument for lsoda

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
    
    ## Get statistics from simulation using a periodic spline interpolant.
    ## The interpolant is fit to the final year of simulation time. 
    ind  <- which(sol$day >= maxtimes-365)
    SeroNegfun <- splinefun(sol[ind,c('day','FracSusc')], method = "periodic")

    ## Store statistics in parmat
    parmat$MinSusfrac[i] <- min(SeroNegfun(0:365))
    parmat$MaxVfrac[i] <- 1 - min(SeroNegfun(0:365))
    parmat$MaxSusfrac[i] <- max(SeroNegfun(0:365))
    parmat$MinVfrac[i] <- 1 - max(SeroNegfun(0:365))
    parmat$AvgSusfrac[i] <- mean(SeroNegfun(0:365))
    parmat$AvgVfrac[i] <- 1 - mean(SeroNegfun(0:365))

    ## Write the data by adding a
    ## the data from row i of parmat
    ## to the data file [filedatname]
    write.table(parmat[i,], file = filedatname,
                append = TRUE, row.names = FALSE, col.names = FALSE)
    
    ## Print and return index to help with debugging
    print(paste("npar: ", i, " / ", npars, ';  Time: ',Sys.time() - starttime, sep=""))
    return(list(parmat[i,], sol))
    
}## Close the function mcl.fun


## Simulate the system. If mclapply throws an error, try
## using the lapply code that is commented out.
starttime <- Sys.time()
sollist <- mclapply(X=1:npars, FUN=mcl.fun, mc.cores=detectCores() - 2)
## sollist <- lapply(X=1, FUN = mcl.fun)
print(paste('Total Time: ', Sys.time() - starttime))

## Read in the simulated data, reorder the rows, and rewrite to file
parmat_filled <- read.table(file = filedatname, header = TRUE)
parmat_filled_ordered <- parmat_filled[order(parmat_filled$npar),]
write.table(parmat_filled_ordered, file = filedatname,
            append = FALSE, row.names = FALSE, col.names = TRUE)


## Choose line colors from palette
cols <- c('orange', 'brown')
cols <- colorRampPalette(cols)(3)

## Load in parmat, specify vectors for d, tv, and tb values
parmat <- read.table(file = paste('Data/', datname,sep = ''), header = T)
parmat <- parmat[order(parmat$npar),]
dvals <- unique(parmat$d)
tvvals <- unique(parmat$tv)
tbvals <- unique(parmat$tb)

## Load in functions that specify differentials
source('Tools/Functions.r', local = FALSE)

fn <- 'Figure_1.eps' ## Filename
setEPS()
postscript(file = fn, height = 3.5, width = 7.08661)
par(mfrow = c(1,2), mar = c(2.5,1.2,1.5,0), oma = c(0, 2, 0, 1))
for(i in c(2,1)){
    wi <- parmat$tb==tbvals[[i]]
    NPeak <- unique(as.numeric(parmat[wi,'NPeak']))
    plot(NA, xlim = c(0, 365), ylim = c(0,1), lwd = 2, xaxt = 'n', yaxt = 'n')
    d <- unique(parmat[wi,'d'])
    tb <- tbvals[[i]]
    polygon(x = c(0,tb,tb,0), y = c(-100,-100,1.03,1.03),
            col = 'gray90', border = NA)
    for(j in 1:length(tvvals)){
        sol <- sollist[[which(wi)[j]]][[2]]
        sol$day <- sol$day%%365
        sol <- sol[-nrow(sol),]
        lines(V/NPop~day, sol, col = cols[j],type = 'l', lwd = 2)
        lines(NPop/NPeak~day,sol, lwd = 0.5, lty = 1, col = 'black')
    }
    xticks <- seq(0,365, by = 30.4)
    xlabs <- c(0,NA,NA,3,NA,NA,6,NA,NA,9,NA,NA,12)
    axis(side = 1, labels = xlabs, at = xticks, padj = -0.7)
    mtext(side = 1, text = 'Time (Month)', line = 1.5)
    yticks <- seq(0,1,by = 0.2)
    ylabs <- yticks
    axis(side = 2, labels = if(i==2){ylabs}else{FALSE}, at = yticks, padj = 0.5)
    if(i==2){
        mtext(side = 2, text = 'Seroprevalence', line = 2.2)
        text(x = 30.4*7, y = 0.90, labels = 'Scaled Population Size', cex = 0.5)
        arrows(x0 = 30.4*7, x1 = 30.4*7, y0 = 0.925, y1 = 0.98, lwd = 1, length = 0.05)
    }
    if(i==1){
        mtext(side = 3, text = '2 Month Birthing Season', line = 0.5)
        text(x = 30.4*7, y = 0.975, labels = 'Scaled Population Size', cex = 0.5)
        arrows(x0 = 30.4*7, x1 = 30.4*6., y0 = 0.955, y1 = 0.885, lwd = 1, length = 0.05)
    }else{
        mtext(side = 3, text = '12 Month Birthing Season', line = 0.5)
    }
    arrows(x0 = tvvals, y0 = 0.1, y1 = 0, col = cols, lwd = 2,
           length = 0.1)
}
dev.off()

