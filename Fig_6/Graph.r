## This script produces a 2x2 paned figure from the Gillespie simulation
## data.

rm(list = ls(all = TRUE))

## Set seed so that bootstrap sampling is reproducible
set.seed(1)

## We use lsoda later on to graph the stable population cycle
require(deSolve) 

## The year that annual vaccination begins. This is the value
## of VaccStartTime in the Gillespie simulation. Set to 2 years
## by default. 
VaccStartYear <- 2 

## Retrieve the parameter matrix corresponding to Gillespie simulation
datname = 'Mastomys'
parname = paste('Gillespie_Simulation/Data/',datname, '/ParMat', sep='')
parmat = read.table(file = parname)
names(parmat) <- c('NPar', 'b', 'd', 'Bp', 'Nv', 'tv', 'gamv', 'gamp', 'tb',
                   'T', 'NPeak')
parmat$Rp <- round(with(parmat, Bp/(d + gamp)))

## Retrieve matrix of extinction times
TExtName <- paste('Gillespie_Simulation/Data/', datname, '/TExtMat', sep='')
TExtMat = read.table(TExtName)

## Specify parameters used in the simulations
tvvals <- unique(parmat$tv)
ntvvals <- length(tvvals)
NvVals <- unique(parmat$Nv)
Rpval <- unique(parmat$Rp)
gamp <- unique(parmat$gamp)
tbval <- unique(parmat$tb)
NPeakval <- unique(parmat$NPeak)

## This script will output a 1x2 figure for each
## row of gridvals, defined here: 
gridvals <- expand.grid(NvVal = NvVals, Rpval = Rpval)
NValues <- nrow(gridvals)

## Years to get reduction / elimination statistics
focYears <- c(1,2,3,4,5) 
nfocYears <- length(focYears)

## Columns to retain from simulation data 
keep <- c('time', 'S', 'Sv', 'V', 'Ip', 'P', 'N')

## Set up arrays that will store the bootstrapped simulation
## results. The arrays ExtBoot and RedBoot, defined below, store
## the bootstrapped extinction probability, and yearly reduction in
## mean pathogen number, for each bootstrapped sample. 
ntrials <- dim(TExtMat)[2]
ExtBoot <- array(data = NA, dim = c(NValues, ntvvals, nfocYears, ntrials))
RedBoot <- array(data = NA, dim = c(NValues, ntvvals, nfocYears, ntrials))

## Loop through each row of gridvals, and produce a 2x2 figure
## given that row's parameters.
for(jj in 1:NValues){

    NvVal <- gridvals$NvVal[jj]

    ## Loop through different timing values, tv
    for(tvi in 1:ntvvals){

        tvval <- tvvals[tvi]
        parrow <- which(parmat$Nv==NvVal & parmat$gamp==gamp & parmat$tv==tvval &
                        parmat$Rp==Rpval)
        wi <- parrow
        NPar <- parmat$NPar[parrow]

        ## Load in the simulation data
        simdat <- read.table(file = paste('Gillespie_Simulation/Data/', datname,
                                          '/Par_', NPar,sep=''), header = TRUE)
        simdat <- simdat[,keep]

        ## Find places in simdat where simulation restarted because of pathogen
        ## eradication or because the maximum simulation time had been reached.
        ## These are places where time column is equal to 0.
        breakpoints <- which(simdat$time==0)
        breakpoints <- c(breakpoints,nrow(simdat))
        nbreakpoints <- length(breakpoints)
        nsims       <- length(breakpoints) - 1
        
        ## These matrices will store the mean pathogen load before (meansimbef)
        ## and after (meansimaft) vaccination has begun. 
        meansimaft <- matrix(, nrow = nsims, ncol = nfocYears)
        meansimbef <- matrix(, nrow = nsims, ncol = 1)
        
        ## For each simulation break, we'll 
        ## get mean abundance of pathogen-infected hosts prior to vaccination.
        ## Start assessment 1 year before vaccination begins
        minfoctime.path <- (VaccStartYear - 2)*365 + tvval 
        ## End assessment immediately before the first pulse vaccination 
        maxfoctime.path <- (VaccStartYear)*365 + tvval 

        ## We'll get the mean abundance of pathogen-infected hosts after
        ## vaccination as well (below), and use these quantities to
        ## calculate the fractional reduction in pathogen load. 
        for(breaki in 1:nsims){

            ## Filter out that part of simdat corresponding to a single simulation
            simdati <- simdat[breakpoints[breaki]:(breakpoints[breaki+1]-1),]
            
            ## Note that we remove any NA values. As a result, simulations that
            ## end before the vaccination has begun are included in this data. This
            ## is pathogen elimination due to stochasticity. We filter these out
            ## in the bootstrapping stage (below).
            simapprox <- approxfun(simdati$time, simdati$Ip)
            meansimbef[breaki,1] <- mean(simapprox(seq(minfoctime.path, maxfoctime.path, by = 1)))

            
            ## Loop through the focal years for which assessment statistics are needed
            for(foci in 1:length(focYears)){
                
                focYear <- focYears[foci]
                
                minfoctime.vacc <- (focYear-1)*365 + tvval ## start assessment at vaccination time tv
                maxfoctime.vacc <- (focYear)*365 + tvval ## end assessment at tv + 365
                
                ##Interpolate and average. 
                maxt <- max(simdati$time)
                simapprox <- approxfun(simdati$time, simdati$Ip)
                meansimaft[breaki,foci] <- mean(simapprox(seq(minfoctime.vacc, maxfoctime.vacc, by = 1)))
            } ##End loop through foci
        } ## End loop through breaki
        
        ## Count extinctions that occured by the end of year focYear
        for(foci in 1:length(focYears)){

            ## We sample from the simulations with replacement. For each sampled set, we calculate 
            ## 1) the average fractional reduction in the mean number of pathogen-infected hosts, and
            ## 2) the probability of elimination status over the set.
            ## We generate a 95\% confidence interval of vaccination outcomes by repeating the
            ## bootstrap procedure. In the line below, we limit our focus to simulations that
            ## have the pathogen present immediately before vaccination begins. 
            valid.sims <- TExtMat[parrow,] >= (tvval + VaccStartYear*365+ tvval)
            total.valid.sims <- sum(valid.sims)
            for(booti in 1:length(valid.sims)){
                valid.boot.set <- which(valid.sims) 
                boot.samp <- sample(valid.boot.set, replace = TRUE, size = total.valid.sims)

                ## Select bootstrap sample
                ExtBoot[jj,tvi,foci,booti] = sum(TExtMat[parrow,boot.samp] <=
                                           (tvval + focYears[foci]*365))/total.valid.sims

                reddat <- 1 - meansimaft[boot.samp, foci]/meansimbef[boot.samp,1]
                RedBoot[jj,tvi, foci, booti] <- mean(reddat, na.rm = TRUE)
            }
        }
        
        print(paste('jj: ', jj,'; tvi:',tvi, '; tv: ', tvval,  '; nsims:', nsims), quote = FALSE)
    }## Close tvi loop
    print(paste("************************** ", jj, " *************************", sep = ''),
          quote = FALSE)
}




## Calculate a 100*(1 - alpha/2)% confidence interval for the values in ExtBoot and RedBoot.  
alpha = 0.025
array.ext = apply(ExtBoot, MARGIN = c(1,2,3),
               FUN = function(x){quantile(x, probs = c(alpha, 0.5, 1 - alpha), na.rm = TRUE)})
array.red = apply(RedBoot, MARGIN = c(1,2,3),
               FUN = function(x){quantile(x, probs = c(alpha, 0.5, 1 - alpha), na.rm = TRUE)})


## Make the figure

## Get population dynamics to plot total population size
source('Gillespie_Simulation/Tools/Functions.r',local=TRUE)
y0 <- c(10000,0,0,0,0)
out <- data.frame(lsoda(y = y0, times = seq(0,10*365,by=1),
                        func = rhs.freq.full, parms = parmat[1,],
                        maxsteps = 100000))
names(out) = c('time', 'S','Sv','Ip','V','P')
out$N = with(out, S + Sv + Ip + V + P)

cx <- 0.6
foc.years = rev(c(1,3,5)) ## Focal years to include in the plot

fn <- 'Figure_6.eps'
setEPS()
postscript(file = fn, height = 5.55, width = 5)
par(mfrow = c(2,2), omi = c(0,0,0,0), mai = c(0.5, 0.4, 0.3, 0.05))
for(jj in 1:NValues){
    ## Proportional Reduction
    plot(NA, xlim = c(0,365), ylim = c(0,1.2), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    polygon(x = c(0,120,120,0), y =c(-100,-100,1.2*1.03,1.2*1.03), col = 'gray90',
            border = NA)
    for(foc.year in foc.years){
        polygon(x = c(tvvals, rev(tvvals)),y = c(array.red[1,jj,,foc.year],
                                                 rev(array.red[3,jj,,foc.year])),
                col = foc.year + 1, border = NA)
        lines(tvvals, array.red[2,jj,,foc.year], lty = 3)
    }
    wi.lastyear <-  max(out$time) - 365 < out$time  & out$time < max(out$time)
    matlines(out$time[wi.lastyear]%%365, out$N[wi.lastyear]/NPeakval)
    legend(x = 'topright',legend = foc.years, col = foc.years + 1, pch = 19,
           title = as.expression(bquote(underline('Year'))),
           cex = cx, bg = 'white')
    nvval <- gridvals$NvVal[jj]
    legend(x = -40, y = 1.25,legend = as.expression(bquote(.(c(500,1000)[jj])~' Vaccinations')),
           cex = cx, bg = 'white')
    axis(side = 1, labels = FALSE, at = seq(0,365, by = 365/12) )
    text(x = -10 + seq(0,365, length.out = 13), -0.125, labels = month.abb[c(6:12,1:6)],
         srt = 65, pos = 1, xpd = TRUE, cex = cx)
    axis(side = 2, labels = TRUE, at = seq(0,1,by = 0.2), cex.axis = cx, padj = 1)
    mtext(side = 3, text = 'Reduction in Pathogen Abundance', line = 0.25, cex = cx)
    mtext(side = 1, text = 'Time of Vaccination', line = 1.8, cex = cx)


    ## Eradication
    plot(NA, xlim = c(0,365), ylim = c(0,1.2), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    polygon(x = c(0,120,120,0), y =c(-100,-100,1.2*1.03,1.2*1.03), col = 'gray90',
            border = NA)
    for(foc.year in foc.years){
        polygon(x = c(tvvals, rev(tvvals)),y = c(array.ext[1,jj,,foc.year], rev(array.ext[3,jj,,foc.year])),
                col = foc.year + 1, border = NA)
        lines(tvvals, array.ext[2,jj,,foc.year], lty = 3)
    }
    wi.lastyear <-  max(out$time) - 365 < out$time  & out$time < max(out$time)
    matlines(out$time[wi.lastyear]%%365, out$N[wi.lastyear]/NPeakval)
    legend(x = 'topright',legend = foc.years, col = foc.years + 1, pch = 19,
           title = as.expression(bquote(underline('Year'))),
           cex = cx, bg = 'white')
    legend(x = -40, y = 1.25,legend = as.expression(bquote(.(c(500,1000)[jj])~' Vaccinations')),
           cex = cx, bg = 'white')
    mtext(side = 3, text = 'Probability of Elimination', line = 0.25, cex = cx)
    axis(side = 2, labels = TRUE, at = seq(0,1,by = 0.2), cex.axis = cx, padj = 1)
    axis(side = 1, labels = FALSE, at = seq(0,365, by = 365/12) )
    text(x = -10 + seq(0,365, length.out = 13), -0.125, labels = month.abb[c(6:12,1:6)],
         srt = 65, pos = 1, xpd = TRUE, cex = cx)
    mtext(side = 1, text = 'Time of Vaccination', line = 1.8, cex = cx)
} ## Close loop through gridvals
dev.off()

