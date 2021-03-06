## Check Gillespie simulation data against ODE simulation. Outputs
## figures to Checkdata folder (in Data folder)

rm(list = ls(all = TRUE))
SimName = "Mastomys"
parmat = read.table(file = paste("Data/",SimName,"/ParMat", sep=''), header = F)
names(parmat) = c('Par','b','d','Bp','Nv','tv','gamv','gamp','tb','T', 'NPeak')

require(deSolve)
source('Tools/Functions.r', local = TRUE)

times = seq(0,10000)
y0 = c(429, 0, 40, 0, 909)
names(y0) = c('S','Sv','Ip','V','P')

endtime <- 365*10
timeseq  <- seq(0,endtime, by = 0.01)

i = 14

Par = parmat$Par[i]
parms = parms = parmat[i,]
VaccPer        <- parmat$T[i]
    NVacc          <- parmat$Nv[i]

    tv <- parmat$tv[i]
    vacctimes <- seq(tv, 6*365,by = VaccPer)
   
    print(paste("NPar:", i-1))
    filename = paste("Data/", SimName, "/Par_", Par, sep = '')
    dat = read.table(file = filename, header = T, sep = " ")

    out <- data.frame(lsoda(y = y0, times = seq(0, endtime, by = 1), 
                            func = rhs.freq.full, parms = parms,
                            events=list(func = vaccinate.full, time=vacctimes),
                            maxsteps = 100000))



##########################################
## WITH VACCINATION
fn  = paste("Data/CheckData/")#Figure folder
dir.create(fn, showWarnings = FALSE)
##Big Picture
    pdf(file = paste(fn, "big_",Par,".pdf",sep = ""), height = 4, width = 6)
    par(mai = c(1,1,0.1,0.1))
    plot(NA, xlim = c(0,endtime), ylim = c(0,1500), type = 'l', lwd = 1)
    these <- seq(1,nrow(dat), by = 100)
    points(S~time, dat[these,], col = 'black', cex = 0.25, pch = 1)
    points(Sv~time, dat[these,], col = 'green', cex = 0.25, pch = 1)
    points(Ip~time, dat[these,], col = 'red', cex = 0.25, pch = 1)
    points(V~time, dat[these,], col = 'lightblue', cex = 0.25, pch = 1)
    points(P~time, dat[these,], col = 'pink', cex = 0.25, pch = 1)

    lines(S~time, out, col = 'gray', lty = 1, lwd = 1)
lines(Sv~time, out, col = 'darkgreen', lty = 1, lwd = 1)
    lines(Ip~time, out, col = 'darkred', lty = 1, lwd = 1)
    lines(V~time, out, col = 'darkblue', lty = 1, lwd = 1)
    lines(P~time, out, col = 458, lty = 1, lwd = 1)
    
    legend = c('S','Sv','Ip','V','P')
    legend(x = "topright", legend = legend, col = c('black', 'green', 'red', 'blue', 'pink'), lwd = 2)
    
    abline(v = vacctimes, lty = 3)

    dev.off()

#########################

    #Small Picture
    pdf(file = paste(fn, "small_",Par,".pdf",sep = ""), height = 4, width = 6)
    par(mai = c(1,1,0.1,0.1))
    plot(NA, xlim = c(0,endtime), ylim = c(0,200), type = 'l', lwd = 1)
    these <- seq(1,nrow(dat), by = 100)
    points(S~time, dat[these,], col = 'black', cex = 0.25, pch = 1)
    points(Sv~time, dat[these,], col = 'green', cex = 0.25, pch = 1)
    points(Ip~time, dat[these,], col = 'red', cex = 0.25, pch = 1)
    points(V~time, dat[these,], col = 'blue', cex = 0.25, pch = 1)
    points(P~time, dat[these,], col = 'pink', cex = 0.25, pch = 1)

    lines(S~time, out, col = 'gray', lty = 1, lwd = 1)
lines(Sv~time, out, col = 'darkgreen', lty = 1, lwd = 1)
    lines(Ip~time, out, col = 'darkred', lty = 1, lwd = 1)
    lines(V~time, out, col = 'darkblue', lty = 1, lwd = 1)
    lines(P~time, out, col = 458, lty = 1, lwd = 1)


    legend = c('S','Sv','Ip','V','P')
    legend(x = "topright", legend = legend, col = c('black', 'green', 'red', 'blue', 'pink'), lwd = 2)
    
    abline(v = vacctimes, lty = 3)

    dev.off()

#########################

##########################################

## Simulation without vaccination
filename = paste("Data/Mastomys_Base", "/Par_", 0, sep = '')
dat.base = read.table(file = filename, header = T, sep = " ")
parms$Nv <- 0
out.base <- data.frame(lsoda(y = y0, times = seq(0, endtime, by = 1), 
                             func = rhs.freq.full, parms = parms,
                             events=list(func = vaccinate.full, time=vacctimes),
                             maxsteps = 100000))


## WITHOUT VACCINATION
fn  = paste("Data/CheckData/")#Figure folder
##Big Picture
    pdf(file = paste(fn, "big_",Par,"_Base.pdf",sep = ""), height = 4, width = 6)
    par(mai = c(1,1,0.1,0.1))
    plot(NA, xlim = c(0,endtime), ylim = c(0,1500), type = 'l', lwd = 1)
    
    these <- seq(1,nrow(dat.base), by = 100)
    points(S~time, dat.base[these,], col = 'black', cex = 0.25, pch = 1)
    points(Sv~time, dat.base[these,], col = 'green', cex = 0.25, pch = 1)
    points(Ip~time, dat.base[these,], col = 'red', cex = 0.25, pch = 1)
    points(V~time, dat.base[these,], col = 'blue', cex = 0.25, pch = 1)
    points(P~time, dat.base[these,], col = 'pink', cex = 0.25, pch = 1)

    lines(S~time, out.base, col = 'gray', lty = 1, lwd = 1)
lines(Sv~time, out.base, col = 'darkgreen', lty = 1, lwd = 1)
    lines(Ip~time, out.base, col = 'darkred', lty = 1, lwd = 1)
    lines(V~time, out.base, col = 'darkblue', lty = 1, lwd = 1)
    lines(P~time, out.base, col = 458, lty = 1, lwd = 1)


    legend = c('S','Sv','Ip','V','P')
    legend(x = "topright", legend = legend, col = c('black', 'green', 'red', 'blue', 'pink'), lwd = 2)
    
    abline(v = vacctimes, lty = 3)

    dev.off()

#########################

    #Small Picture
    pdf(file = paste(fn, "small_",Par,"_Base.pdf",sep = ""), height = 4, width = 6)
    par(mai = c(1,1,0.1,0.1))
    plot(NA, xlim = c(0,endtime), ylim = c(0,200), type = 'l', lwd = 1)
    
    these <- seq(1,nrow(dat.base), by = 100)
    points(S~time, dat.base[these,], col = 'black', cex = 0.25, pch = 1)
    points(Sv~time, dat.base[these,], col = 'green', cex = 0.25, pch = 1)
    points(Ip~time, dat.base[these,], col = 'red', cex = 0.25, pch = 1)
    points(V~time, dat.base[these,], col = 'blue', cex = 0.25, pch = 1)
    points(P~time, dat.base[these,], col = 'pink', cex = 0.25, pch = 1)

    lines(S~time, out.base, col = 'gray', lty = 1, lwd = 1)
lines(Sv~time, out.base, col = 'green', lty = 1, lwd = 1)
    lines(Ip~time, out.base, col = 'red', lty = 1, lwd = 1)
    lines(V~time, out.base, col = 'blue', lty = 1, lwd = 1)
    lines(P~time, out.base, col = 'pink', lty = 1, lwd = 1)


    legend = c('S','Sv','Ip','V','P')
    legend(x = "topright", legend = legend, col = c('black', 'green', 'red', 'blue', 'pink'), lwd = 2)
    
    abline(v = vacctimes, lty = 3)

    dev.off()

#########################
