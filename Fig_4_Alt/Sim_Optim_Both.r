rm(list = ls(all = TRUE))

## Simulation time: 316 cpu minutes (when tb's length.out
## set to 25. Figure generated with length.out = 365).

SimName = "Fig_4_alt"
parmat = expand.grid(NPar = 0, b = NA,
       d = c(1/(225), 1/365),
      	             rho = c(0.5),
                     gamv = 1/14, tb = seq(1,365, length.out = 25),
       trans = c('freq', 'dens'),
       T = 365, NPeak = 1000, tvstar = 0,
                     downtv80 = 0, uptv80 = 0,
                     downtv90 = 0, uptv90 = 0,
                     downtv95 = 0, uptv95 = 0,
                     downtv975 = 0, uptv975 = 0,
		     vavgstar = 0, vnull = 0)

trans = 'freq' ## Transmission type

NPars = nrow(parmat)
parmat$NPar = 1:NPars
parmat$b     <- with(parmat,
                     NPeak*d/(exp(d*(T-tb))*(exp(d*tb)-1)/(exp(d*T)-1)))
parmat$Nv <- with(parmat, NPeak*rho)

require(deSolve)
require(parallel)
source('Tools/Functions.r', local = TRUE)

nyears   <- 100
maxtimes <- 365*nyears
times    <- seq(0,maxtimes, by = 1)
witimes <- times >= (nyears-1)*365
yeartimes <- c(times[witimes]%%365)
yeartimes[length(yeartimes)] <- 365
tvvals <- matrix(c(0.1,seq(1,364.9,length.out = 25)), ncol = 1)


filename = paste("Data/", SimName, sep = '')
write.table(parmat[1,][-1,], file = filename, append = F)

mcl.fun = function(i){
    cost.fxn <- function(tv,parms){
        TVal  <- parmat$T[i]
        NVacc <- parms$Nv
        vacctimes <- seq(tv, maxtimes,by = TVal)
        y0 = c(1000, 0, 0)
        names(y0) = c('S','Sv','V')
        out  <- data.frame(lsoda(y = y0, times = times,
                                 func = rhs.jv, parms = parms,
                                 events=list(func = vaccinate.jv, time=vacctimes),
                                 maxsteps = 100000))
###Return average serorpevalence
        ## Change to density dependence assessment

if(trans=='dens'){
        vals    <- with(out[witimes,], V)
	vals[1] <- vals[length(vals)]
	vapprox = splinefun(x = yeartimes, vals, method = "periodic")
        VAvg <- 1/365*(integrate(vapprox, lower = 0, upper = 365, stop.on.error = FALSE)$value)
	vals    <- with(out[witimes,], S + Sv + V)
	vals[1] <- vals[length(vals)]
	NPopapprox = splinefun(x = yeartimes, vals, method = "periodic")
        NAvg <- 1/365*(integrate(NPopapprox, lower = 0, upper = 365, stop.on.error = FALSE)$value)
        return(VAvg/NAvg)
}else{
    vals    <- with(out[witimes,], V/(S + Sv + V))
    vals[1] <- vals[length(vals)]
    vapprox = splinefun(x = yeartimes, vals, method = "periodic")
    VAvg <- 1/365*(integrate(vapprox, lower = 0, upper = 365, stop.on.error = FALSE)$value)
    return(VAvg)
}


    }###End cost.fxn

    vnullvals <- apply(X = tvvals, MARGIN = 1, FUN = cost.fxn, parms = parmat[i,])

    wi.interval <- c(max(1, which.max(vnullvals) - 1), min(which.max(vnullvals) + 1, length(tvvals)))
    interval <- tvvals[wi.interval]

    ## Find optimal time of vaccination and corresponding average seroprevalence
    uroot <- optimize(f=cost.fxn, interval = interval, parmat[i,], maximum=TRUE)
    tvstar <- uroot$maximum ###Best tv value
    vavgstar <- uroot$objective ###Maximum value

    ## Find and store 90% and 95% optimal strategy
    wi.down <- which(tvvals <= tvstar)
    wi.up <- which(tvvals >= tvstar)

    if(length(wi.down) > 1){
        f.down.approx <- approxfun((vnullvals/vavgstar)[wi.down], tvvals[wi.down])
        ## 90% optimal
        if(!is.na(f.down.approx(0.9))){
            parmat$downtv90[i] <- f.down.approx(0.9)
        }else{
            parmat$downtv90[i] <- 0
        }
        ## 95% optimal
        if(!is.na(f.down.approx(0.95))){
            parmat$downtv95[i] <- f.down.approx(0.95)
        }else{
            parmat$downtv95[i] <- 0
        }
    }else{
        parmat$downtv90[i] <- 0
        parmat$downtv95[i] <- 0
    }

    if(length(wi.up) > 1){
        f.up.approx <- approxfun((vnullvals/vavgstar)[wi.up], tvvals[wi.up])
        if(!is.na(f.up.approx(0.9))){
        parmat$uptv90[i] <- f.up.approx(0.9)
        }else{
            parmat$uptv90[i] <- 365
        }
        if(!is.na(f.up.approx(0.95))){
            parmat$uptv95[i] <- f.up.approx(0.95)
        }else{
            parmat$uptv95[i] <- 365
        }
    }else{
        parmat$uptv90[i] <- 0
        parmat$uptv95[i] <- 0
    }


    parmat$tvstar[i] <- tvstar
    parmat$vavgstar[i] <- vavgstar

    vnullfun = splinefun(x = c(tvvals,365.1), c(vnullvals,vnullvals[1]))
    parmat$vnull[i]  <- 1/365*integrate(vnullfun, lower = 0, upper = 365)$value

    ## write data
    write.table(parmat[i,], file = filename, append = T, col.names = F, row.names = F)
    print(paste("NPar:", i, ' / ', NPars, ' ; Time: ', Sys.time() - starttime),
          quote = FALSE)
    return(i)
}
starttime <- Sys.time()
mcl.out <- mclapply(X = 1:NPars, FUN = mcl.fun, mc.cores = detectCores()-2)

refparmat <- read.table(file = filename, header = T)
refparmat <- refparmat[order(refparmat$NPar),]
write.table(refparmat, file = filename, append = F, row.names = F)
parmat <- read.table(file = filename, header = T)

