rm(list = ls(all=TRUE))
datname = 'Fig_2'
filedatname <- paste("Data/",datname,sep="")
parmat = read.table(file = filedatname, header = T)

XValName = 'tv'
YValName = 'VdivNPop' ## 'AvgVfrac'
FValName = 'tb'
ZValName = 'd'

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
for(type in c('.eps')){
    fn <- paste('Figure_2',type, sep='')

    setEPS()
    postscript(file = fn, height = 5, width = 5)
layout(matrix(c(1,2,
                3,4), nrow = 2, byrow = T), widths = c(0.4,0.4,0.0),
       heights = c(0.25,0.25,0.25,0.25))
    par(oma = c(0,2,0,0.5), mar = c(3, 1.5, 1.5, 0))
for(i.f in 1:nFVals){
    FVal = FVals[i.f]
    plot.mat = matrix(nrow = nXVals, ncol = nZVals)
    for(i.z in 1:nZVals){
        wi <- parmat[,ZValName] == ZVals[i.z] & parmat[,FValName] == FVals[i.f]
        plot.mat[,i.z] = parmat[wi,YValName]
    }
    plot.mat = apply(plot.mat, 2, function(x){x/max(x)}) 
    matplot(NA, type = 'l', xlim = c(0,365), ylim = c(0.90,1),
            xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    polygon(x = c(0.0, FVal, FVal, 0.0), y = c(-2,-2,1.03,1.03), col = 'gray90',
            lty = 3, lwd = 0.5, border = NA)
    matlines(XVals, plot.mat, lwd = 1.5, lty = 1, col = cols)
    mm <- apply(X = plot.mat, MARGIN = 2, FUN = which.max)
    mmind <- cbind(mm,c(1,2,3,4))
##    matpoints(XVals[mm], plot.mat[mmind], pch = 1)
#    abline(v = FVals[i.f], lwd = 0.25, lty = 3)

###Axes
    yticks = seq(0,1,by = 0.1)
    if(i.f %in% c(1,3)){
        axis(side = 2, labels = yticks,
             at = yticks, padj = 0.5, cex.axis = 0.7)
mtext(side = 2, text = 'Scaled Seroprevalence', line = 2.25, outer = FALSE, cex = 0.7)

    }else{
        axis(side = 2, labels = FALSE, at = yticks,
             padj = 0.5)
    }


    xticks = round(seq(0,365, by = 30.4))
    labticks = c(0, NA, NA, 3, NA, NA, 6, NA, NA, 9, NA, NA, 12)
###Axis labels
    if(i.f %in% c(1,2)){
        axis(side = 1, labels = FALSE, at = xticks)
    }
    if(i.f %in% c(1,2,3,4)){
        mtext(side = 1, text = 'Month of Vaccination', line = 1.4, cex = 0.7)
        axis(side = 1, labels = labticks, at = xticks, padj = -1.25, cex.axis = 0.7)
    }

    if(i.f==1){
        legendtext = rev(c('250 Days', '2.5 Yrs', '5 Yrs', '10 Yrs'))
        legend(x = 'topright', legend = legendtext, col = rev(cols),
               lwd = 2, bty = 'n', cex = 0.7,
               title = as.expression(bquote(underline('Lifespan'))))

    }

    mtext(side = 3, text =  paste('Birthing Season: ',
                                  c('1 Month', '3 Months', '6 Months','9 Months')[i.f],sep=''),
          cex = 0.7)
}
dev.off()
}## End loop through figure types
