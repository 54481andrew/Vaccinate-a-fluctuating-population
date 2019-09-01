
rm(list = ls(all = T))
require(deSolve)

datname = 'Fig_3'
filedatname <- paste("Data/",datname,sep="")
parmat = read.table(file = filedatname, header = T)

XValName = 'tv'
YValName = 'AvgVfrac'
ZValName = 'rho'

XVals = unique(parmat[,XValName])
FVals = c(1,2)
ZVals = unique(parmat[,ZValName])

nXVals = length(XVals)
nFVals = length(FVals)
nZVals = length(ZVals)

##cols <- c('darkgreen', 'lightblue')
cols <- c('blue', 'orange')
cols <- colorRampPalette(cols)(3)
source('Tools/Functions.r', local = T)
tb <- unique(parmat$tb)

### Get simulation data for total population
### size.
    pars = parmat[1,]
    times <- seq(0,365*10, by = 1)
    yinit <- c(100,0,0)
    out       <- data.frame(lsoda(y = yinit, times = times,
                       hmax = 0.1, func = rhs.jv,
                       parms = pars))
    names(out) <- c('day', 'S', 'Sv', 'V')
    out$NPop <- with(out, S + Sv + V)
    NPeak <- pars$NPeak
    out <- subset(out, subset = day >= 9*365)
    out <- out[-nrow(out), ]

###

### Graph both file types
for(type in c('.eps')){
    fn <- paste('Figure_3',type, sep='')

    setEPS()
    postscript(file = fn, height = 3.14961, width = 3.14961)
par(mar = c(2.75,1.2,1,0), oma = c(0, 2, 0, 1))
    YValName <- 'AvgVfrac'
    plot.mat = matrix(nrow = nXVals, ncol = nZVals)
    for(i.z in 1:nZVals){
        wi <- parmat[,ZValName] == ZVals[i.z]
            plot.mat[,i.z] = parmat[wi,YValName]
    }
    matplot(NA, type = 'l', xlim = c(0,365), ylim = c(0,1),
            xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    polygon(x = c(0.0, tb, tb, 0.0), y = c(-2,-2,1.03,1.03), col = 'gray90',
            border = NA)
##    segments(x0 = 0, y0 = -2, y1 = 3, lty = 3, lwd = 0.5)
##    segments(x0 = tb, y0 = -2, y1 = 3, lty = 3, lwd = 0.5)
    matlines(XVals, plot.mat[,3:1], lwd = 2, lty = 1, col = rev(cols))
###Axes
    yticks = seq(0,1,by = 0.1)
        axis(side = 2, labels = yticks,
             at = yticks, padj = 0.5, cex.axis = 0.75)
    axis(side = 1, labels = seq(0,12, by = 3), at = seq(0,365, by = 30.4*3), padj = -0.6 )
###Axis labels
    mtext(side = 1, text = 'Month of Vaccination', line = 1.85, cex = 1)
        legendtext = rev(c('250', '500', '750'))
        legend(x = 30.4*7, y = 1, legend = legendtext, col = rev(cols),
               lwd = 2, bty = 'n', cex = 0.7,
               title = as.expression(bquote(underline('Vaccinations'))))
    ## Add line that shows population size, scaled to max size
    pars = parmat[1,]
    times <- seq(0,365*10, by = 1)
    yinit <- c(100,0,0)
    out       <- data.frame(lsoda(y = yinit, times = times,
                       hmax = 0.1, func = rhs.jv,
                       parms = pars))
    names(out) <- c('day', 'S', 'Sv', 'V')
    out$NPop <- with(out, S + Sv + V)
    NPeak <- pars$NPeak
    out <- subset(out, subset = day >= 9*365)
    out <- out[-nrow(out), ]
    lines(NPop/NPeak~(day%%365),out, lwd = 0.5, lty = 1, col = 'black')
    f <- approx((out$NPop/NPeak)[out$day%%365 > tb], out$day[out$day%%365 > tb], xout = c(0.25,0.5,0.75))
    ts <- f$y
   ### matpoints(ts%%365, c(0.25, 0.5, 0.75), pch = 1)
mtext(side = 2, text = 'Average Seroprevalence', line = 1.05, outer = TRUE, cex = 1)
dev.off()
} ## End loop through figure types