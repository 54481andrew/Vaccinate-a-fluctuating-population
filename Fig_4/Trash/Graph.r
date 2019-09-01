## This script reads the data written by Simulate.r, and
## graphs the fractional reduction in mean pathogen abundance.
## The plot is a 2x2 panel figure.

rm(list = ls(all=TRUE))

## Read in data
datname = 'Fig_4' ## Name of datafile to be read in 
filedatname <- paste("Data/",datname,sep="")
parmat = read.table(file = filedatname, header = T)
parmat$f.red = with(parmat, 1 - Ipv/Ip0) ## calculate fractional reduction

## Specify the structure of the graph. 
XValName = 'tv' ## Varied along x-axis of plots
YValName = 'f.red' ## Varied along y-axis of plots
FVal1Name = 'type' ## Varied across different panes in the horizontal direction
FVal2Name = 'pmu' ## Varied across different panes in the vertical direction
ZValName = 'gamp' ## Different lines in same pane

XVals = unique(parmat[,XValName])
FVal1s = unique(parmat[,FValName])
ZVals = unique(parmat[,ZValName])

nXVals = length(XVals)
nFVals = length(FVals)
nZVals = length(ZVals)

initcols = c('blue', 'orange')
cols = colorRampPalette(c('blue', 'orange'))(length(ZVals))

## source('Tools/Functions.r', local = T)
parmat$d <- round( parmat$d , 5)
parmat$gamp <- round( parmat$gamp , 5)
parmat$gamv <- round( parmat$gamv , 5)

## Parameters fixed over entire figure
dvals <- unique(parmat$d)
gamvvals <- unique(parmat$gamv)
gampvals <- unique(parmat$gamp)
Rpvals <- unique(parmat$Rp)
tbvals <- unique(parmat$tb)
rhovals <- unique(parmat$rho)

gridvals <- expand.grid(d = dvals, Rp = Rpvals, tb = tbvals, gamp = gampvals,
                        gamv = gamvvals, rho = rhovals)

for(ii in 1:nrow(gridvals)){
    
    dval <- gridvals[ii,'d']
    Rpval <- gridvals[ii,'Rp']
    tbval <- gridvals[ii,'tb']
    gamvval <- gridvals[ii,'gamv']
    gampval <- gridvals[ii,'gamp']
    rhoval <- gridvals[ii,'rho']


    
    ## Graph both file types
    for(type in c('.png')){
    
        wi.fig <- parmat$d == dval & parmat$Rp == Rpval &
            parmat$gamv == gamvval & parmat$gamp == gampval &
            parmat$rho== rhoval ## parameters to keep in the fig
        
        dn <- round(1/dval/365,1)
        gpn <- round(1/gampval)
        gvn <- round(1/gamvval)

        dir.create(paste('Figures/LS_',dn, '_rho_' , rhoval , sep = ''), showWarnings = FALSE)
        
        fn <- paste( 'Figures/', paste('LS_',dn, '_rho_' , rhoval , sep = ''), '/',
                    'Rp_', Rpval,
                    '_gamp_', gpn, '_gamv_', gvn, 
                    '.png', sep = '' )
        
        ##    setEPS()
        ##    postscript(file = fn, height = 5, width = 5)
        png(file = fn, height = 3, width = 6, res = 400, units = 'in')
        layout(matrix(c(1,2), nrow = 1, byrow = T), widths = c(0.4,0.4),
               heights = c(0.25,0.25))
        par(oma = c(0,2,0,0.5), mar = c(3, 1.5, 1.5, 0))
        for(i.f in 1:nFVals){
            FVal = FVals[i.f]
            plot.mat = matrix(nrow = nXVals, ncol = nZVals)
            for(i.z in 1:nZVals){
                wi <- parmat[,ZValName] == ZVals[i.z] & parmat[,FValName] == FVals[i.f] & wi.fig
                plot.mat[,i.z] = parmat[wi,YValName]
            }
            wi.maxs <- apply(plot.mat, 2, function(x){XVals[which.max(x)]})
            max.vals <-apply(plot.mat, 2, max) 
            matplot(NA, type = 'l', xlim = c(0,365), ylim = c(0.0,1),
                    xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
            polygon(x = c(0.0, tbval, tbval, 0.0), y = c(-2,-2,1.03,1.03), col = 'gray90',
                    lty = 3, lwd = 0.5, border = NA)
            matlines(XVals, plot.mat, lwd = 1.5, lty = 1, col = cols)
            points(wi.maxs, max.vals, pch = 19, cex = 0.5, col = cols)
            ## Axes
            if(i.f==1){
                yticks = seq(0,1,by = 0.1)
                axis(side = 2, labels = yticks,
                     at = yticks, padj = 0.5, cex.axis = 0.7)
                mtext(side = 2, text = 'Reduction in Mean Infected', line = 2.25,
                      outer = FALSE, cex = 0.7)
            }else{
                axis(side = 2, labels = FALSE,
                     at = yticks, padj = 0.5, cex.axis = 0.7)
            }
            
            xticks = round(seq(0,365, by = 30.4))
            labticks = c(0, NA, NA, 3, NA, NA, 6, NA, NA, 9, NA, NA, 12)
            ## Axis labels
            
            axis(side = 1, labels = FALSE, at = xticks)
            mtext(side = 1, text = 'Month of Vaccination', line = 1.4, cex = 0.7)
            axis(side = 1, labels = labticks, at = xticks, padj = -1.25, cex.axis = 0.7)
            
            
            mtext(side = 3, text = ifelse(parmat$type[wi] == 'f' , 'Frequency-Dependent',
                                          'Density-Dependent'))
        
            if(i.f==1){
                legendtext = as.expression(sapply(1:length(ZVals),
                                                  function(i){bquote('p'[mu]*' = '*.(ZVals[i]))}))
                legend(x = 'topright', legend = legendtext, col = rev(cols),
                       lwd = 2, bty = 'n', cex = 0.7,
                       title = as.expression(bquote(underline('Virulence'))))
            }
        }
        dev.off()
    } ## End loop through figure types
} ## Loop through gridvals (paramters that are fixed in the figure)