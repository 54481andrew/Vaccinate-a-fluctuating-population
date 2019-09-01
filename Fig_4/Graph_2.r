## This script reads the data written by Simulate.r, and
## graphs the fractional reduction in mean pathogen abundance.
## The plot is a 1x2 panel figure.

rm(list = ls(all=TRUE))

## Read in data
datname = 'Fig_4' ## Name of datafile to be read in 
filedatname <- paste("Data/",datname,sep="")
parmat = read.table(file = filedatname, header = T)
parmat$f.red = with(parmat, 1 - Ipv/Ip0) ## calculate fractional reduction

## Specify the structure of the graph. 
XValName = 'tv' ## Varied along x-axis of plots
YValName = 'f.red' ## Varied along y-axis of plots
FVal1Name = 'type' ## Different panes
ZVal1Name = 'Rp' ## Different lines in same pane
ZVal2Name = 'gamp' ## Different lines in same pane

## Round some parameters to help with certain steps below
parmat$d <- round( parmat$d , 5)
parmat$gamp <- round( parmat$gamp , 5)
parmat$gamv <- round( parmat$gamv , 5)

## Extract values to plot
XVals = unique(parmat[,XValName])
FVals1 = unique(parmat[,FVal1Name])
ZVals1 = unique(parmat[,ZVal1Name])
ZVals2 = unique(parmat[,ZVal2Name])
nXVals = length(XVals)
nFVals1 = length(FVals1)
nZVals1 = length(ZVals1)
nZVals2 = length(ZVals2)

## Specify colors used in the plot
col1 = 'blue'
col2 = 'orange'

## Parameters fixed over entire figure
dvals <- unique(parmat$d)
gamvvals <- unique(parmat$gamv)
Rpvals <- unique(parmat$Rp)
tbvals <- unique(parmat$tb)
Nvvals <- unique(parmat$Nv)
pmuvals <- unique(parmat$pmu)

gridvals <- expand.grid(d = dvals, tb = tbvals,
                        gamv = gamvvals, Nv = Nvvals,
                        pmu = pmuvals)

for(ii in 1:nrow(gridvals)){

    dval <- gridvals[ii,'d']
    tbval <- gridvals[ii,'tb']
    gamvval <- gridvals[ii,'gamv']
    Nvval <- gridvals[ii,'Nv']
    pmuval <- gridvals[ii,'pmu']
    
    ## Which rows of parmat correspond to the fixed parameters in this figure
    wi.fig <- parmat$d == dval &
        parmat$gamv == gamvval & parmat$Nv == Nvval &
        parmat$pmu == pmuval

    ## Write an informative filename
    dn <- round(1/dval/365,1) ## Get lifespan in years
    gvn <- round(1/gamvval) ## Get gamv in 
    dir.create(paste('Figures_1/', paste('LS_',dn, sep = ''), sep = ''), showWarnings = FALSE) 
    fn <- paste( 'Figures_1/', paste('LS_',dn, sep = ''), '/',
                'pmu_', pmuval,'gamv_', gvn, 'Nv_', Nvval, 
                '.eps', sep = '' )
    
    setEPS()
    postscript(file = fn, height = 3, width = 6)
    layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE), widths = c(0.5,0.5),
           heights = c(0.5,0.5))
    par(oma = c(2.5,2.1,0,0.5), mar = c(0.5, 1, 0.5, 0))

    ## Loop through each figure pane and plot the correct x-y relationship
    for(i.f1 in 1:nFVals1){
            FVal1 = FVals1[i.f1]
            plot.mat = matrix(nrow = nXVals, ncol = nZVals1*nZVals2)
            for(i.z1 in 1:nZVals1){
                for(i.z2 in 1:nZVals2){
                    wi <- parmat[,ZVal1Name] == ZVals1[i.z1] &
                        parmat[,ZVal2Name] == ZVals2[i.z2] &
                        parmat[,FVal1Name] == FVals1[i.f1] &
                        wi.fig
                    plot.mat[,2*(i.z1-1) + i.z2] = parmat[wi,YValName]
                }
            }
            
            
            wi.maxs <- apply(plot.mat, 2, function(x){XVals[which.max(x)]})
            max.vals <-apply(plot.mat, 2, max)
            matplot(NA, type = 'l', xlim = c(0,365), ylim = c(0.0,1.),
                    xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
            polygon(x = c(0.0, tbval, tbval, 0.0), y = c(-2,-2,0.9,0.9), col = 'gray90',
                    lty = 3, lwd = 0.5, border = NA)
            polygon(x = c(0.0, tbval, tbval, 0.0), y = c(0.6,0.6,1.01,1.01), col = 'white',
                    lty = 3, lwd = 0.5, border = NA)
            matlines(XVals, plot.mat, lwd = 1.5, pch = c(1,1,2,2), lty = c(1,2,1,2),
                     type = 'l', col = c(col1, col1, col2, col2))
            points(wi.maxs, max.vals, pch = 19, cex = 0.5, col = 'black')
            
            ## y axis
            if(i.f1==1){
                yticks = seq(0,1,by = 0.1)
                axis(side = 2, labels = yticks,
                     at = yticks, padj = 0.5, cex.axis = 0.7)
                mtext(side = 2, text = 'Reduction in Mean Infected', line = 2.25,
                      outer = FALSE, cex = 0.7)
            }else{
                axis(side = 2, labels = FALSE,
                     at = yticks, padj = 0.5, cex.axis = 0.7)
            }

            ## x axis
        xticks = round(seq(0,365, by = 30.4))
        labticks = c(0, NA, NA, 3, NA, NA, 6, NA, NA, 9, NA, NA, 12)
        
        axis(side = 1, labels = FALSE, at = xticks)        
        text(x = 180, y = 0.95, labels = ifelse(parmat$type[wi] == 'f' , 'Frequency-Dependent',
                                                 'Density-Dependent'))

            mtext(side = 1, text = 'Month of Vaccination', line = 1.4,
                  cex = 0.7)
            axis(side = 1, labels = labticks, at = xticks,
                 padj = -1.25, cex.axis = 0.7)
            


            if(i.f1==1){
                
                legendtext = as.expression(sapply(1:length(ZVals1),
                                                  function(i){bquote('R'['0,p']*' = '*.(ZVals1[i]))}))
                legend(x = -5, y = 0.85, legend = legendtext, col = c(col1,col2),
                       lwd = 2, bty = 'n', cex = 0.7,
                       title = as.expression(bquote(underline('Pathogen Transmission'))))
                
                
                legendtext = as.expression(sapply(1:length(ZVals2),
                                                  function(i){bquote(gamma[p]*' = '*.(round(1/ZVals2[i])))}))
                legend(x = 0.55*365, y = 0.85, legend = legendtext, col = 'black',
                       lwd = 2, bty = 'n', cex = 0.7, lty = c(1,2),
                       title = as.expression(bquote(underline('Pathogen Recovery'))))
                
            }
    } ## Loop through different figure panes
    dev.off()
} ## Loop through gridvals (paramters that are fixed in the figure)
    
