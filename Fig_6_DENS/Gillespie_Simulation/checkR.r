## Check Gillespie simulation data against ODE simulation. Outputs
## figures to Checkdata folder (in Data folder)

rm(list = ls(all = TRUE))
SimName = "Mastomys"
parmat = read.table(file = paste("Data/",SimName,"/ParMat", sep=''), header = F)
names(parmat) = c('Par','b','d','Bp','Nv','tv','gamv','gamp','tb','T', 'NPeak')
parmat$pmu = 0

require(deSolve)
source('Tools/Functions.r', local = TRUE)

times = seq(0,10000)
y0 = c(497, 0, 29, 0, 852)
names(y0) = c('S','Sv','Ip','V','P')

VaccStartTime <- 2*365
endtime <- 365*10
timeseq  <- seq(0,endtime, by = 0.01)

i = 1

Par = parmat$Par[i]
parms = parms = parmat[i,]
VaccPer        <- parmat$T[i]
    NVacc          <- parmat$Nv[i]

    tv <- parmat$tv[i]
    vacctimes <- seq(tv + VaccStartTime, 6*365,by = VaccPer)
   
    print(paste("NPar:", i-1))
    filename = paste("Data/", SimName, "/Par_", Par, sep = '')
    dat = read.table(file = filename, header = T, sep = " ")

    out <- data.frame(lsoda(y = y0, times = seq(0, endtime, by = 1), 
                            func = rhs.dens.full, parms = parms,
                            events=list(func = vaccinate.full, time=vacctimes),
                            maxsteps = 100000))



##########################################
## WITH VACCINATION
fn  = paste("Data/CheckData/")#Figure folder
dir.create(fn, showWarnings = FALSE)
##Big Picture
    pdf(file = paste(fn, "overview_",Par,".pdf",sep = ""), height = 4, width = 6)
    par(mai = c(1,1,0.1,0.1))
    plot(NA, xlim = c(0,endtime), ylim = c(0,1750), type = 'l', lwd = 1)
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
    pdf(file = paste(fn, "closeup_",Par,".pdf",sep = ""), height = 4, width = 6)
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
