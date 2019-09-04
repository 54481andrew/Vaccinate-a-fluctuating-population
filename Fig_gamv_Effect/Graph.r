rm(list = ls(all=TRUE))
datname = 'Fig_S'
filedatname <- paste("Data/",datname,sep="")
parmat = read.table(file = filedatname, header = T)

XValName = 'tv'
YValName = 'AvgVfrac'
FValName = 'd'
ZValName = 'gamv'

XVals = unique(parmat[,XValName])
FVals = unique(parmat[,FValName])
ZVals = unique(parmat[,ZValName])
ZVals = ZVals[]


nXVals = length(XVals)
nFVals = length(FVals)
nZVals = length(ZVals)

initcols = c('blue', 'orange')
cols = colorRampPalette(c('blue', 'orange'))(length(ZVals))

source('Tools/Functions.r', local = T)

### Graph both file types
    fn <- 'Figure_S.eps'

    setEPS()
    postscript(file = fn, height = 5, width = 5)
    par(mar = c(3, 3, 1.5, 0.5))
i.f = 1
    FVal = FVals[i.f]
    plot.mat = matrix(nrow = nXVals, ncol = nZVals)
    for(i.z in 1:nZVals){
        wi <- parmat[,ZValName] == ZVals[i.z] & parmat[,FValName] == FVals[i.f]
        plot.mat[,i.z] = parmat[wi,YValName]
    }
    matplot(NA, type = 'l', xlim = c(0,365), ylim = c(0.0,1),
            xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    polygon(x = c(0.0, 120, 120, 0.0), y = c(-2,-2,1.03,1.03), col = 'gray90',
            lty = 3, lwd = 0.5, border = NA)
    matlines(XVals, plot.mat, lwd = 1.5, lty = 1, col = cols)
    mm <- apply(X = plot.mat, MARGIN = 2, FUN = which.max)
    mmind <- cbind(mm,c(1,2,3,4))

###Axes
    yticks = seq(0,1,by = 0.1)
    if(i.f %in% c(1,3)){
        axis(side = 2, labels = yticks,
             at = yticks, padj = 0.5, cex.axis = 0.7)
mtext(side = 2, text = 'Average Seroprevalence', line = 2.25, outer = FALSE, cex = 0.7)

    }else{
        axis(side = 2, labels = FALSE, at = yticks,
             padj = 0.5)
    }

    xticks = round(seq(0,365, by = 30.4))
    labticks = c(0, NA, NA, 3, NA, NA, 6, NA, NA, 9, NA, NA, 12)
    ## Axis labels
    if(i.f %in% c(1,2)){
        axis(side = 1, labels = FALSE, at = xticks)
    }
    if(i.f %in% c(1,2,3,4)){
        mtext(side = 1, text = 'Month of Vaccination', line = 1.4, cex = 0.7)
        axis(side = 1, labels = labticks, at = xticks, padj = -1.25, cex.axis = 0.7)
    }

    if(i.f==1){
        legendtext = (c('0 Days', '7 days', '14 days', '30 days', '60 days'))
        legend(x = 0.06*365, y = 0.45, legend = legendtext, col = cols,
               lwd = 2, bty = 'n', cex = 0.7,
               title = as.expression(bquote(underline('Immunity Lag'))))

    }

    mtext(side = 3, text =  paste('Birthing Season: ',
                                  '4 Months',sep=''),
          cex = 0.7)

dev.off()

