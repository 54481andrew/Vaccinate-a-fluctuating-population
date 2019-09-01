## This script reads the data written by Simulate.r, and
## graphs the optimal timing of vaccination. The plot is a 2x2 paned
## figure.

rm(list = ls(all = TRUE))

## Read in data
datname = "Fig_3" ## Name of datafile to be read in 
filedatname <- paste("Data/",datname,sep="")
parmat = read.table(file = filedatname, header = T)

## Specify the structure of the graph. 
XValName = 'tb' ## Varied along x-axis
YValName = 'tvstarsc' ## Varied along y-axis
FValName = 'd' ## Varied across different panes 

## Extract values to plot
XVals = unique(parmat[,XValName])
FVals = unique(parmat[,FValName])
nXVals = length(XVals)
nFVals = length(FVals)

## Specify colors used in the plot
require(RColorBrewer)
cols = brewer.pal(9, 'Blues')
cols = cols[c(3,5)]

## File to which the figure is saved
fn <- 'Figure_3.eps'

setEPS()
postscript(fn, height = 5, width = 5)
layout(matrix(c(1,2,3,4), nrow = 2, byrow = T), widths = c(0.5,0.5),
       heights = c(0.5,0.5))
par(oma = c(2.5,2.1,0,0.5), mar = c(0.5, 1, 0.5, 0))

## Loop through 4 different FVals, each of which corresponds to a single
## pane. 
for(i.f in 1:nFVals){
    FVal = FVals[i.f]

    ## Each column of plot.mat will contain a vaccination window.
    ## These are used to plot the optimal timing, as well as
    ## a polygon around the optimal timing that indicates
    ## the vaccination window. 
    plot.mat = matrix(nrow = nXVals, ncol = 5) 
    plot.mat.loess = plot.mat
    wi <- parmat[,FValName] == FVals[i.f]
    plot.mat[,1] = parmat[wi,'downtv90']
    plot.mat[,2] = parmat[wi,'downtv95']
    plot.mat[,3] = parmat[wi,'tvstar']
    plot.mat[,4] = parmat[wi,'uptv95']
    plot.mat[,5] = parmat[wi,'uptv90']
    
    matplot(NA, type = 'l', xlim = c(0,365), ylim = c(-30,400),
            xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    
    ## 90% Border
    polygon(x = c(XVals, rev(XVals)), y = c(plot.mat[,1], rev(plot.mat[,5])),
            col = cols[1], border = NA)
    ## 95% Border
    polygon(x = c(XVals, rev(XVals)), y = c(plot.mat[,2], rev(plot.mat[,4])),
            col = cols[2], border = NA)

    these <- seq(1,length(XVals)-1,by = 3)
    lines(XVals[these], plot.mat[these,3], lwd = 1.5, lty = 1, col = 'black', pch = 8, cex = 1,
          type = 'b')

    ## Axes
    yticks = seq(0,365, length.out = 13)
    yticklabs = c(0,NA,2,NA,4,NA,6,NA,8,NA,10,NA,12)
    if(i.f %in% c(1,3)){
        axis(side = 2, labels = yticklabs,
             at = yticks, padj = 0.5, cex.axis = 1)
    }else{
        axis(side = 2, labels = FALSE, at = yticks,
             padj = 0.5)
    }

    xticks = seq(0,365, by = 30)
    xlabticks = c(0, NA,2,NA,4,NA,6,NA,8,NA,10,NA,12)

    ## Axis labels
    if(i.f %in% c(1,2)){
        axis(side = 1, labels = FALSE, at = xticks)
    }
    if(i.f %in% c(3,4)){
        mtext(side = 1, text = 'Birthing Duration (Months)', line = 1.75, cex = 1)
        axis(side = 1, labels = xlabticks, at = xticks, cex.axis = 1,
             padj = -0.75)
    }

    xlabpos <- 180
    ylabpos <- 395
	    text(x = xlabpos, y = ylabpos, labels =
	    	   bquote(underline(.(round(1/(365*FVal),2))~' Year Lifespan')), cex = 1)

    ## Descriptive legend in panel 4
    if(i.f==4){
        legendtext = c('Optimal Timing', '> 90% Optimal', '> 95% Optimal')
        legend(x = 'bottomright', legend = legendtext, pch = c(8,22,22), lty = c(2,NA,NA),
               bg = 'white',
               col = c('black', 'black', 'black'), pt.bg = c(NA, cols[1], cols[2]),
               pt.cex = c(1, 1.75, 1.75), lwd = c(2, NA, NA), pt.lwd = c(1,0.5,0.5),
               cex = 0.8)
    }
} ## End loop through i.f values (panes)

mtext(side = 2, text = 'Optimal Time of Vaccine Delivery (Month)', line = 1.05, outer = TRUE, cex = 1)
dev.off()
