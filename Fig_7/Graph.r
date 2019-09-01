## This script simulates Lassa virus circulating in a
## population of Meles meles

############################################
############# USER INPUT ###################
############################################

rm(list = ls(all=TRUE))
ParMat <- expand.grid(d = c(1/(4*365)),
                      gamv = 1/14,
                      Nv = c(250),
                      tb = 60.8, T = 365, NPeak = 1107,
                      gamp = c(1/(4*365)), pmu = c(1),
		      Rp = c(2),
                      tv = c(30, 120, 180))

## Use the specified peak population size to
## determine the constant birth rate during
## the breeding season. This relationship
## NPeak and b is derived in the supplemental
## Mathematical file.
ParMat$b     <- with(ParMat,
                     NPeak*d/(exp(d*(T-tb))*(exp(d*tb)-1)/(exp(d*T)-1)))
## Calculate Bp from Rp
ParMat$Bp    <- with(ParMat, Rp*(d+gamp))

## Add simulation with no vaccination
ParMat <- rbind(ParMat, c(ParMat[1,]))
ParMat[nrow(ParMat),'Nv'] <- 0

############################################
############## END USER INPUT ##############
############################################

## Load packages for lsoda numerical
## integration. Functions.r for functions that
## return differentials of the ODE.
require(deSolve)
source("Tools/Functions.r", local = T)

## Give each row of parmat a name
npars <- nrow(ParMat)
ParMat$npar <- 1:npars

## Specify times to simulate for. TBurnIn period
## allows pathogen to settle into stable cycles. maxtimes
## is the end time of the simulation.
TBurnIn <- 20*365 ## Pathogen circulates for 20 years
maxtimes <- 40*365 ## Max simulation time set to 40 years

## Loop through each parameter set, numerically integrate
## the ODE system, then store the results to meta.out.
meta.out <- list()
starttime <- Sys.time()
for(i in 1:nrow(ParMat)){
    print(i)
    tvval <- ParMat$tv[i]
    tv.first <- TBurnIn + tvval ### Get first vaccination time

    ## Initialize system with 5 infected hosts, everyone else susceptible
    y0 = c( with(ParMat[i,], NPeak*exp(-d*(365-tb))), 0, 5, 0, 0)
    names(y0) = c('S','Iv','Ip','V','P')
    NVacc          <- ParMat$Nv[i]
    ParMat$Nv[i] <- 0 ## temporarily set to 0 during "no vaccination" phase of simulation
    out1 <- data.frame(lsoda(y = y0, times = seq(0,TBurnIn,by=1),
                             func = rhs.freq.full, parms = ParMat[i,]))
    y0 = out1[nrow(out1), -1]
    y0 = as.numeric(y0)
    ParMat$Nv[i] <- NVacc
    names(y0) = c('S','Iv','Ip','V','P')
    vacctimes <- seq(TBurnIn + tvval, maxtimes, by = 365)
    out2 <- data.frame(lsoda(y = y0, times = seq(TBurnIn, maxtimes, by = 1),
                             func = rhs.freq.full, parms = ParMat[i,],
                             events=list(func = vaccinate.full, time=vacctimes),
                             maxsteps = 100000))
    out = rbind(out1, out2[-1,])
    meta.out[[i]] <- out
    print(paste('Simulation', i, 'complete'), quote = FALSE)
}

## Simulation without vaccination
parms <- ParMat[1,]
parms$Nv <- 0
parms$b     <- with(parms,
                     NPeak*d/(exp(d*(T-tb))*(exp(d*tb)-1)/(exp(d*T)-1)))
parms$Bp    <- 0
y0 = with(parms, c(NPeak*exp(-d*(365-tb)), 0, 5, 0, 0))
names(y0) = c('S','Iv','Ip','V','P')
NVacc          <- parms$Nv[i]
source("Tools/Functions.r", local = T)
out.base <- data.frame(lsoda(y = y0, times = seq(0,maxtimes,by=1),
                        func = rhs.freq.full, parms = parms))
names(out.base) = c('time', 'S','Iv','Ip','V','P')
out.base$N <- with(out.base, S + Iv + Ip + V + P)


############################

## Graph Figure that shows dynamics of each class, for each of the
## 3 vaccination times.

fn <- 'Figure_S4.eps'
setEPS()
    postscript(file = fn, height = 5, width = 5)
par(omi = c(0.5,0.5,0.1,0.1), mai = c(0.1, 0.0, 0.0, 0.0))
layout(matrix(1:3, ncol = 1), respect = FALSE)
for(i in 1:(nrow(ParMat)-1)){
    tvval <- ParMat$tv[i]
    tv.first <- TBurnIn + tvval ### Get first vaccination time
    out <- data.frame(meta.out[[i]])
    tv <- ParMat$tv[i]
    tv.first <- tv + TBurnIn
    matplot(NA, xlim = c(TBurnIn-2*365, TBurnIn + 5*365), ylim = c(0,1200),
            xaxt = 'n', yaxt = 'n',
            xlab = '', ylab = '', )
    for(yi in 1:50){
        time <- yi*365
        polygon(x = c(time, time + ParMat$tb[i], time + ParMat$tb[i], time),
                y = c(-5000, -5000, 1200*1.03, 1200*1.03), col = 'gray90',
                border = NA)
    }
    matlines(out$time, cbind(out$S, out$Ip),
             col = c('blue', 'red'), lty = 1, lwd = 1.5, type = 'l')
    taxis <- seq(0,365*50, by = 3*365/12)
    if(i==3){labs = 1:50 - 10; mtext(side = 1, text = 'Time (Year)', line = 2)}else{labs = F}
    axis(side = 1, labels = labs, at = (1:50)*365, tick = TRUE,
         cex.axis = 0.75, padj =-1)
    axis(side = 2,
    labels = c(0, NA, 500, NA, 1000),
    at = seq(0,1000, by = 250), cex.axis = 0.75, padj = 1)
    mtext(side = 2,
    text = as.expression(bquote('Number of Hosts')),
    line = 2., outer = TRUE)
    arrows(x0 = tv.first + seq(0,4380, by = 365), y0 = 1200, y1 = 1000,
           lwd = 1.5, col = 'orange', length = 0.075)
    matlines(out.base$time, out.base$N, lty = 1, col = 'black')
    if(i==1){
        legend(x = 'topleft', legend = c('Total','S', 'Ip'), col = c('black', 'blue', 'red', 'gray'),
               lwd = 2, lty = 1, cex = 0.75, ncol = 1, bg = 'white')
    }
}
dev.off()


############################

## Graph figure that shows number of pathogen-infected hosts
## on the y-axis, time on the x-axis, for each of the 3 times
## of vaccination. Add a line for the case of no vaccination
## too.

cols <- c('orange', 'red', 'purple', 'gray') ## Line colors
fn <- 'Figure_7.eps'
setEPS()
postscript(file = fn, height = 3, width = 5)
par(omi = c(0,0.0,0,0), mai = c(0.5, 0.55, 0.1, 0.05))
matplot(NA, xlim = c(TBurnIn - 1.3*365, TBurnIn + 3*365),
        ylim = c(0,600), xaxt = 'n', yaxt = 'n',
        xlab = '', ylab = '')
for(yi in 1:50){
    time <- yi*365
    polygon(x = c(time, time + ParMat$tb[i], time + ParMat$tb[i], time),
            y = c(-100, -100, 600*1.03, 600*1.03), col = 'gray90',
            border = NA)
}
for(i in 1:nrow(ParMat)){
    matlines(out$time, meta.out[[i]][,'Ip'],
             col = cols[i], lty = 1, lwd = 1.5, type = 'l')
    taxis <- seq(0,365*50, by = 3*365/12)
    tv.first <- TBurnIn + ParMat$tv[i]
    if(i < nrow(ParMat)){
        arrows(x0 = tv.first + seq(0,2*365, by = 365), y0 = 595.5, y1 = 540,
               lwd = 1.5, col = cols[i], length = 0.1)
    }
}
axis(side = 1, labels = FALSE, at = taxis)
axis(side = 1, labels = 1:50 - 10, at = (1:50)*365, tick = F,
     cex.axis = 0.75, padj =-1)
axis(side = 2, labels = TRUE, at = seq(0,1000, by = 100), cex.axis = 0.75,
     padj = 1)
mtext(side = 1, text = 'Time (Year)', line = 1.5)
mtext(side = 2,
      text = 'Number of Infected Hosts',
      line = 1.5)
tvvals <- unique(ParMat$tv)
legend(x = 'topleft', legend = c(tvvals,'No Intervention'),
       col = cols, lwd = 2, lty = 1,
       cex = 0.75,
       title = as.expression(bquote(underline('Vaccination Day'))), bg = 'white')
dev.off()
