rm(list = ls(all=TRUE))

ParMat <- expand.grid(d = c(1/(225)),
                      gamv = 1/14,
                      rho = c(0.1855),
                      tb = 120, T = 365, b = 16.67,##NPeak = 2696,
                      gamp = c(1/30.4), pmu = c(0.0),
		      Rp = c(2),
                      tv = c(60, 120, 180))
name <- paste('freq_ls',round(1/unique(ParMat$d)), sep = '')

##ParMat$b     <- with(ParMat,
##                     NPeak*d/(exp(d*(T-tb))*(exp(d*tb)-1)/(exp(d*T)-1)))
ParMat$NPeak <- with(ParMat, b*exp(d*(T-tb))*(exp(d*tb) - 1)/(d*(exp(d*T) - 1)))
## CHANGE BP WHEN CHANGING d
ParMat$Bp    <- 0.07467836##0.0712689 ##with(ParMat, Rp*(d+gamp))
ParMat$Nv    <- with(ParMat, NPeak*rho)

## Simulation with no vaccination
ParMat <- rbind(ParMat, c(ParMat[1,]))
ParMat[nrow(ParMat),'Nv'] <- 0

############################################
############## END USER INPUT ##############
############################################

require(deSolve)
require(parallel)



npars <- nrow(ParMat)
ParMat$npar <- 1:npars
TBurnIn <- 10*365 ## Pathogen circulates for 10 years
maxtimes <- 20*365 ## Max simulation time

meta.out <- list()
starttime <- Sys.time()
for(i in 1:nrow(ParMat)){
    print(i)
    tvval <- ParMat$tv[i]
    tv.first <- TBurnIn + tvval ### Get first vaccination time
    y0 = c(500*exp(-0.01*(365-120)), 0, 5, 0, 0)
    names(y0) = c('S','Iv','Ip','V','P')
    NVacc          <- ParMat$Nv[i]
    ParMat$Nv[i] <- 0 #temporarily set to 0
    source("Tools/Functions.r", local = T)
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
    print('done')
}

### Simulation without vaccination
parms <- ParMat[1,]
parms$Nv <- 0
parms$b     <- with(parms,
                     NPeak*d/(exp(d*(T-tb))*(exp(d*tb)-1)/(exp(d*T)-1)))
parms$Bp    <- 0
y0 = with(parms, c(NPeak*exp(-d*(365-tb)), 0, 0, 0, 0))
names(y0) = c('S','Iv','Ip','V','P')
NVacc          <- parms$Nv[i]
source("Tools/Functions.r", local = T)
out.base <- data.frame(lsoda(y = y0, times = seq(0,maxtimes,by=1),
                        func = rhs.freq.full, parms = parms))
names(out.base) = c('time', 'S','Iv','Ip','V','P')
out.base$N <- with(out.base, S + Iv + Ip + V + P)



#################
## COMBINE TIME##
#################

### Graph eps file type
for(type in c('.png')){
    fn <- paste('Figures/mFigure_S3_', name,type, sep='')
    ##    setEPS()
##    postscript(file = fn, height = 5, width = 5)
    png(file = fn, height = 5, width = 5, units = 'in', res = 400)
    par(omi = c(0.5,0.5,0.1,0.1), mai = c(0.1, 0.0, 0.0, 0.0))
layout(matrix(1:3, ncol = 1), respect = FALSE)
for(i in 1:(nrow(ParMat)-1)){
    tvval <- ParMat$tv[i]
    tv.first <- TBurnIn + tvval ### Get first vaccination time
    out <- data.frame(meta.out[[i]])
    tv <- ParMat$tv[i]
    tv.first <- tv + TBurnIn
    plot(NA, , xlim = c(TBurnIn-2*365, TBurnIn + 5*365), ylim = c(0,5000), xaxt = 'n', yaxt = 'n',
         xlab = '', ylab = '')
    for(yi in 1:20){
        time <- yi*365
        polygon(x = c(time, time + ParMat$tb[i], time + ParMat$tb[i], time),
                y = c(0, 0, 36050, 36050), col = 'gray90',
                border = NA)
    }
    matlines(out$time, cbind(out$S, out$Ip, out$Ip + out$P),
                             col = c('blue', 'red', 'darkgray'), lty = 1, lwd = 1.5, type = 'l')
    taxis <- seq(0,365*20, by = 3*365/12)
    if(i==3){labs = 1:20 - 10; mtext(side = 1, text = 'Time (Year)', line = 2)}else{labs = F}
    ##axis(side = 1, labels = FALSE, at = taxis)
    axis(side = 1, labels = labs, at = (1:20)*365, tick = TRUE,
         cex.axis = 0.75, padj =-1)
    axis(side = 2,
    labels = c(0, NA, 10, NA, 20, NA, 30),
    at = seq(0,30000, by = 5000), cex.axis = 0.75, padj = 1)
    mtext(side = 2,
    text = 'Number of Hosts',
    line = 2., outer = TRUE)
##as.expression(bquote('Density ( x10'^'3'~'per km'^'2'~')')),
          arrows(x0 = tv.first + seq(0,4380, by = 365), y0 = 5000, y1 = 4500,
           lwd = 1.5, col = 'orange', length = 0.075)
    matlines(out.base$time, out.base$N, lty = 1, col = 'black')
    if(i==1){
        legend(x = 'topleft', legend = c('Total','S', 'Ip', 'Ip + P'), col = c('black', 'blue', 'red', 'gray'),
               lwd = 2, lty = 1, cex = 0.75, ncol = 2, bg = 'white')
    }
}
dev.off()
} ## Loop through figure types


################
## METAFIGURE ##
################

cols <- c('orange', 'red', 'purple', 'gray')
### Graph both file types
for(type in '.png'){##c('.eps')){
    fn <- paste('Figures/mFigure_5_',name, type, sep='')
##    setEPS()
##    postscript(file = fn, height = 3, width = 5)
    png(file = fn, height = 5, width = 5, units = 'in', res = 400)
    par(omi = c(0,0.0,0,0), mai = c(0.5, 0.55, 0.1, 0.05))
matplot(NA, xlim = c(TBurnIn - 1*365, TBurnIn + 3*365),
        ylim = c(0,260), xaxt = 'n', yaxt = 'n',
        xlab = '', ylab = '')
for(yi in 1:20){
    time <- yi*365
    polygon(x = c(time, time + ParMat$tb[i], time + ParMat$tb[i], time),
            y = c(-100, -100, 200*1.03, 200*1.03), col = 'gray90',
            border = NA)
##    abline(v = c(time, time + ParMat$tb[i]), lty = 3, lwd = 0.25)
}
    for(i in 1:nrow(ParMat)){
    matlines(out$time, meta.out[[i]][,'Ip'],
             col = cols[i], lty = 1, lwd = 1.5, type = 'l')
    taxis <- seq(0,365*20, by = 3*365/12)
    tv.first <- TBurnIn + ParMat$tv[i]
    if(i < nrow(ParMat)){
        arrows(x0 = tv.first + seq(0,2*365, by = 365), y0 = 250, y1 = 230,
               lwd = 1.5, col = cols[i], length = 0.1)
        }
}
    axis(side = 1, labels = FALSE, at = taxis)
    axis(side = 1, labels = 1:20 - 10, at = (1:20)*365, tick = F,
         cex.axis = 0.75, padj =-1)
    axis(side = 2, labels = TRUE, at = seq(0,500, by = 50), cex.axis = 0.75,
         padj = 1)
mtext(side = 1, text = 'Time (Year)', line = 1.5)
mtext(side = 2,
text = 'Number of Infected Hosts',
line = 1.5)
tvvals <- unique(ParMat$tv)
    legend(x = 'topleft', legend = c(tvvals,'No Intervention'), col = cols, lwd = 2, lty = 1,
           cex = 0.7,
           title = as.expression(bquote(underline('Vaccination Day'))), bg = 'white')
text(x = TBurnIn + 2*365, y = 350,labels = name)
    dev.off()
} ## Loop through figure types
