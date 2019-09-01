
## Graph probability of extinction, fractional reduction
rm(list = ls(all = TRUE))
require(deSolve)
require(parallel)
datname = 'Mastomys'
parname = paste('Gillespie_Simulation/Data/',datname, '/ParMat', sep='')
ParMat = read.table(file = parname)
names(ParMat) <- c('NPar', 'b', 'd', 'Bp', 'Nv', 'tv', 'gamv', 'gamp', 'tb',
                   'T', 'NPeak')
ParMat$Rp <- round(with(ParMat, Bp/(d + gamp)))
ParMat$Rho <- with(ParMat, Nv/NPeak)

datname.base = 'Mastomys_Base'
parname.base = paste('Gillespie_Simulation/Data/',datname.base, '/ParMat', sep='')
ParMat_Base = read.table(file = paste(parname.base, sep=''))
names(ParMat_Base) <- c('NPar', 'b', 'd', 'Bp', 'Nv', 'tv', 'gamv', 'gamp', 'tb',
                   'T', 'NPeak')

TExtName <- paste('Gillespie_Simulation/Data/', datname, '/TExtMat', sep='')
TExtMat = read.table(TExtName)

## Form the plotting matrix for Pr(Erad|Year)
## Loop through Month, Year,
tvvals <- unique(ParMat$tv)
ntvvals <- length(tvvals)

Rpvals <- unique(ParMat$Rp)
Rhovals <- unique(ParMat$Rho)
gridvals <- expand.grid(Rhoval = Rhovals, Rpval = Rpvals)
NValues <- nrow(gridvals)

focYears <- c(1,2,3,4,5)
nfocYears <- length(focYears)
cols = c('blue', 'yellow', 'green', 'red', 'gray', 'black')
ltys = c(3,3,3,1,1,1)
keep <- c('time', 'S', 'Iv', 'V', 'Ip', 'P', 'N')

MeanMatList <- list()
ExtMatList <- list()

for(jj in 1:NValues){
    rhoval <- gridvals$Rhoval[jj]
    Rpval  <- gridvals$Rpval[jj]
    gamp <- unique(ParMat$gamp)
    tbval <- unique(ParMat$tb)
    NPeakval <- unique(ParMat$NPeak)
    
    PlotMeanMat <- matrix(,nrow = ntvvals, ncol = nfocYears)
    PlotMeanMatODE <- matrix(,nrow = ntvvals, ncol = nfocYears)
    PlotExtMat <- matrix(,nrow = ntvvals, ncol = nfocYears)
    
### loop and collect baseline data
    simdat.base <- read.table(file = paste('Gillespie_Simulation/Data/',datname.base,'/Par_',0,sep=''), header = TRUE)
    simdat.base <- simdat.base[,keep]
    meansimbef  <- mean(simdat.base$Ip) ## Get average pathogen abundance from simulation without vaccination

    for(tvi in 1:ntvvals){
    tvval <- tvvals[tvi]
    parrow <- which(ParMat$Rho==rhoval & ParMat$gamp==gamp & ParMat$tv==tvval &
                    ParMat$Rp==Rpval)
    wi <- parrow
    NPar <- ParMat$NPar[parrow]
    simdat <- read.table(file = paste('Gillespie_Simulation/Data/', datname, '/Par_', NPar,sep=''), header = TRUE)
    simdat <- simdat[,keep]
    breakpoints <- which(simdat$time==0)
    breakpoints <- c(breakpoints,nrow(simdat))
    nbreakpoints <- length(breakpoints)
    tv.first <- unique(ParMat$TVaccStart[wi]) + tvval ### Get first vaccination time

    nsims       <- length(breakpoints) - 1
    meansimaft <- matrix(, nrow = nsims, ncol = nfocYears)
    extmat <- matrix(0, nrow = nsims, ncol = nfocYears)
    
    ### loop and collect sim data
    for(breaki in 1:nsims){

        simdati <- simdat[breakpoints[breaki]:(breakpoints[breaki+1]-1),]
               
        for(foci in 1:length(focYears)){
            focYear <- focYears[foci]
            minfoctime <- (focYear-1)*365
            maxfoctime <- focYear*365
            wi.focYear <- minfoctime < simdati$time & simdati$time < maxfoctime

            ##Interpolate and average
            maxt <- max(simdati$time)
            simapprox <- approxfun(simdati$time, simdati$Ip)
            f <- function(x){
                ifelse(x > maxt, 0, simapprox(x))
            }
            meansimaft[breaki,foci] <- mean(f(seq(minfoctime, maxfoctime, by = 1)))
        }##foci
    }##breaki
    
    ## Count extinctions that occured by the end of year focYear,
    ## starting at time tv.first
    for(foci in 1:length(focYears)){
        PlotExtMat[tvi,foci] = sum(TExtMat[parrow,] <
                                   focYears[foci]*365)/ncol(TExtMat)
    }

    reddat <- 1 - meansimaft/meansimbef

     PlotMeanMat[tvi,] <- colMeans(reddat, na.rm = T)
    print(paste('jj: ', jj,'; tvi:',tvi, '; tv: ', tvval,  '; nsims:', nsims))
}## Close tvi loop

    MeanMatList[[jj]] <- PlotMeanMat
    ExtMatList[[jj]] <- PlotExtMat
    print(paste("************************** ", jj, " *************************", sep = ''))
}


### Make the figure
tc <- 70
    source('Gillespie_Simulation/Tools/Functions.r',local=TRUE)
    ### Get population dynamics 
    y0 <- c(10000,0,0,0,0)
    out <- data.frame(lsoda(y = y0, times = seq(0,10*365,by=1),
                             func = rhs.freq.full, parms = ParMat_Base[1,],
                             maxsteps = 100000))
    names(out) = c('time', 'S','Iv','Ip','V','P')
    out$N = with(out, S + Iv + Ip + V + P)

cx <- 0.6
### Graph both file types
for(type in c('.eps')){
    fn <- paste('Figure_6',type, sep='')
    setEPS()
    postscript(file = fn, height = 5.55, width = 5)
par(mfrow = c(2,2), omi = c(0,0,0,0), mai = c(0.5, 0.4, 0.3, 0.05))
for(jj in 1:NValues){
    PlotMeanMat <- MeanMatList[[jj]]
    PlotExtMat <- ExtMatList[[jj]]
    ## Proportional Reduction
    plot(NA, xlim = c(0,365), ylim = c(0,1.2), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
polygon(x = c(0,120,120,0), y =c(-100,-100,1.2*1.03,1.2*1.03), col = 'gray90',
        border = NA)
matlines(tvvals, PlotMeanMat, pch = 19, cex = cx, type = 'p',
        lty = 1, lwd = 2)
wi.lastyear <-  max(out$time) - 365 < out$time  & out$time < max(out$time)
matlines(out$time[wi.lastyear]%%365, out$N[wi.lastyear]/NPeakval)
abline(v = seq(0,365, by = 30.4), lwd = 0.5, lty = 3)
abline(h = seq(0,1, by = 0.1), lwd = 0.5, lty = 3)
legend(x = 'topright',legend = 1:5, col = 1:5, pch = 19,
       title = as.expression(bquote(underline('Year'))),
       cex = cx, bg = 'white')
    rhoval <- gridvals$Rhoval[jj]
    legend(x = -40, y = 1.25,legend = as.expression(bquote(.(c(500,1000)[jj])~' Vaccinations')),
           cex = cx, bg = 'white')
    axis(side = 1, labels = FALSE, at = seq(0,365, by = 365/12) )
text(x = -10 + seq(0,365, length.out = 13), -0.125, labels = month.abb[c(6:12,1:6)],
     srt = 65, pos = 1, xpd = TRUE, cex = cx)
    axis(side = 2, labels = TRUE, at = seq(0,1,by = 0.2), cex.axis = cx, padj = 1)
mtext(side = 3, text = 'Reduction in Average Abundance', line = 0.25, cex = cx)
    mtext(side = 1, text = 'Time of Vaccination', line = 1.8, cex = cx)
    ## Eradication
    plot(NA, xlim = c(0,365), ylim = c(0,1.2), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    polygon(x = c(0,120,120,0), y =c(-100,-100,1.2*1.03,1.2*1.03), col = 'gray90',
            border = NA)
    matlines(tvvals, PlotExtMat, pch = 19, cex = cx, type = 'p',
             lty = 1, lwd = 2)
    wi.lastyear <-  max(out$time) - 365 < out$time  & out$time < max(out$time)
    matlines(out$time[wi.lastyear]%%365, out$N[wi.lastyear]/NPeakval)
    abline(v = seq(0,365, by = 30.4), lwd = 0.5, lty = 3)
    abline(h = seq(0,1, by = 0.1), lwd = 0.5, lty = 3)
    legend(x = 'topright',legend = 1:5, col = 1:5, pch = 19,
           title = as.expression(bquote(underline('Year'))),
           cex = cx, bg = 'white')
    legend(x = -40, y = 1.25,legend = as.expression(bquote(.(c(500,1000)[jj])~' Vaccinations')),
           cex = cx, bg = 'white')
    ##mtext(side = 1, text = 'Day of year', line = 2.5)
    mtext(side = 3, text = 'Probability of Elimination', line = 0.25, cex = cx)
    axis(side = 2, labels = TRUE, at = seq(0,1,by = 0.2), cex.axis = cx, padj = 1)
    axis(side = 1, labels = FALSE, at = seq(0,365, by = 365/12) )
    text(x = -10 + seq(0,365, length.out = 13), -0.125, labels = month.abb[c(6:12,1:6)],
     srt = 65, pos = 1, xpd = TRUE, cex = cx)
    mtext(side = 1, text = 'Time of Vaccination', line = 1.8, cex = cx)
} ## Close loop through figure panes
dev.off()
}
