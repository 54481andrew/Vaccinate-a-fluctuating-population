require(parallel)
SimName = "Animals.Fine"

source("make_parmat.r", local = TRUE)

NPars = nrow(parmat)
parmat$NPar = 1:NPars
parmat$b     <- with(parmat,
                     NPeak*d/(exp(d*(T-tb))*(exp(d*tb)-1)/(exp(d*T)-1)))
parmat$Nv <- with(parmat, NPeak*rho)

require(deSolve)
source('Tools/Functions.r', local = TRUE)

nyears   <- 100
maxtimes <- 365*nyears
times    <- seq(0,maxtimes, by = 1)
witimes <- times >= (nyears-1)*365
yeartimes <- c(times[witimes]%%365)
yeartimes[length(yeartimes)] <- 365
tvvals <- matrix(seq(0.1,364.9,length.out = 50), ncol = 1)

filename = paste("Data/", SimName, sep = '')
write.table(parmat[1,][-1,], file = filename, append = F)


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
        
###Return average serorpevalence
	vals    <- with(out[witimes,], V/(S + Sv + V))
	vals[1] <- vals[length(vals)]
	vapprox = splinefun(x = yeartimes, vals, method = "periodic")
        VAvg <- 1/365*(integrate(vapprox, lower = 0, upper = 365, stop.on.error = FALSE)$value)
        return(VAvg)
    }###End cost.fxn

    vnullvals <- apply(X = tvvals, MARGIN = 1, FUN = cost.fxn, parms = parmat[i,])
    tvstar <- tvvals[which.max(vnullvals)]
    
    wi.interval <- c(max(1, which.max(vnullvals) - 1), min(which.max(vnullvals) + 1, length(tvvals)))
    interval <- tvvals[wi.interval]

    uroot <- optimize(f=cost.fxn, interval = interval, parmat[i,], maximum=TRUE)
    tvstar <- uroot$maximum ###Location
    vavgstar <- uroot$objective ###Maximum value

    wi80 <- which(vnullvals/vavgstar > 0.8)
    parmat$downtv80[i] <- min(tvvals[wi80])
    parmat$uptv80[i] <- max(tvvals[wi80])

    wi90 <- which(vnullvals/vavgstar > 0.9)
    parmat$downtv90[i] <- min(tvvals[wi90])
    parmat$uptv90[i] <- max(tvvals[wi90])

    wi95 <- which(vnullvals/vavgstar > 0.95)
    parmat$downtv95[i] <- min(tvvals[wi95])
    parmat$uptv95[i] <- max(tvvals[wi95])

    wi975 <- which(vnullvals/vavgstar > 0.975)
    parmat$downtv975[i] <- min(tvvals[wi975])
    parmat$uptv975[i] <- max(tvvals[wi975])

    
    
    parmat$tvstar[i] <- tvstar
    parmat$vavgstar[i] <- vavgstar

    vnullfun = splinefun(x = c(tvvals,365.1), c(vnullvals,vnullvals[1]))
    parmat$vnull[i]  <- 1/365*integrate(vnullfun, lower = 0, upper = 365)$value

    write.table(parmat[i,], file = filename, append = T, col.names = F, row.names = F)
    print(paste("NPar:", i, ' / ', NPars))
    return(i)
}

metaout <- mclapply(X = 1:NPars, FUN = mcl.fun, mc.cores = detectCores()-2)

refparmat <- read.table(file = filename, header = T)
refparmat <- refparmat[order(refparmat$NPar),]
write.table(refparmat, file = filename, append = F, row.names = F)
parmat <- read.table(file = filename, header = T)

parmat$Window80 <- with(parmat, uptv80 - downtv80)
parmat$Window90 <- with(parmat, uptv90 - downtv90)
parmat$Window95 <- with(parmat, uptv95 - downtv95)
parmat$Window975 <- with(parmat, uptv975 - downtv975)

write.table(parmat, file = filename, append = F, row.names = F)

tab <- subset(parmat, select = c('Animal', 'Window95'))
tab$Month95 <- round(tab$Window95/30.4,1)
tab[,c(1,3)]
