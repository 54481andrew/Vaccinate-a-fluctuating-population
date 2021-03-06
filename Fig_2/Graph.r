## This script reads the data written by Simulate.r, and
## graphs vaccine effectiveness for different timings
## of vaccination. The plot is a 2x2 pane.

rm(list = ls(all=TRUE))

## Read in data
datname = 'Fig_2' ## Name of datafile to be read in 
filedatname <- paste("Data/",datname,sep="")
parmat = read.table(file = filedatname, header = T)

## Specify the structure of the graph. 
XValName = 'tv' ## varied along x-axis
YValName = 'AvgVfrac' ## varied along y-axis
FValName = 'tb' ## fixed within each pane, varied across panes
ZValName = 'd' ## varied across different lines in each pane

## Extract values to plot
XVals = unique(parmat[,XValName])
FVals = unique(parmat[,FValName])
ZVals = unique(parmat[,ZValName])
ZVals = ZVals[]
nXVals = length(XVals)
nFVals = length(FVals)
nZVals = length(ZVals)

## Specify colors used in the plot
initcols = c('blue', 'orange')
cols = colorRampPalette(c('blue', 'orange'))(length(ZVals))

## File to which the figure is saved
fn <- 'Figure_2.eps'

setEPS()
postscript(file = fn, height = 5, width = 5)
layout(matrix(c(1,2,
                3,4), nrow = 2, byrow = T), widths = c(0.4,0.4,0.0),
       heights = c(0.25,0.25,0.25,0.25))
par(oma = c(0,2,0,0.5), mar = c(3, 1.5, 1.5, 0))

## Loop through 4 different FVals, each of which corresponds to a single
## pane. 
for(i.f in 1:nFVals){
    FVal = FVals[i.f]
    plot.mat = matrix(nrow = nXVals, ncol = nZVals)
    for(i.z in 1:nZVals){
        wi <- parmat[,ZValName] == ZVals[i.z] & parmat[,FValName] == FVals[i.f]
        plot.mat[,i.z] = parmat[wi,YValName]
    }
    matplot(NA, type = 'l', xlim = c(0,365), ylim = c(0.0,1),
            xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    polygon(x = c(0.0, FVal, FVal, 0.0), y = c(-2,-2,1.03,1.03), col = 'gray90',
            lty = 3, lwd = 0.5, border = NA)
    matlines(XVals, plot.mat, lwd = 1.5, lty = 1, col = cols)
    mm <- apply(X = plot.mat, MARGIN = 2, FUN = which.max)
    mmind <- cbind(mm,c(1,2,3,4))

    ## Axes
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
        legendtext = rev(c('1 Yr', '2.5 Yrs', '5 Yrs', '10 Yrs'))
        legend(x = 0.06*365, y = 0.45, legend = legendtext, col = rev(cols),
               lwd = 2, bty = 'n', cex = 0.7,
               title = as.expression(bquote(underline('Lifespan'))))
        
    }

    ## Text at top of panes
    mtext(side = 3, text =  paste('Birthing Season: ',
                                  c('1 Month', '3 Months', '6 Months','9 Months')[i.f],sep=''),
          cex = 0.7)
}

dev.off()

