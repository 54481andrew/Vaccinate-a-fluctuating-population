## This script reads the data written by Simulate.r, and
## graphs the fractional reduction in mean pathogen abundance.
## The plot is a 2x2 panel figure.

rm(list = ls(all=TRUE))

## Read in data
datname = 'Fig_4' ## Name of datafile to be read in 
filedatname <- paste("Data/",datname,sep="")
parmat = read.table(file = filedatname, header = T)
parmat$f.red = with(parmat, 1 - Ipv/Ip0) ## calculate fractional reduction

XValName = 'tv' ## Varied along x-axis of plots
YValName = 'f.red' ## Varied along y-axis of plots
FValName = 'type' ## Different panes
ZVal1Name = 'Rp' ## Different lines in same pane
ZVal2Name = 'gamp' ## Different lines in same pane


## source('Tools/Functions.r', local = T)
parmat$d <- round( parmat$d , 5)
parmat$gamp <- round( parmat$gamp , 5)
parmat$gamv <- round( parmat$gamv , 5)

XVals = unique(parmat[,XValName])
FVals = unique(parmat[,FValName])
ZVals1 = unique(parmat[,ZVal1Name])
ZVals2 = unique(parmat[,ZVal2Name])
nXVals = length(XVals)
nFVals = length(FVals)
nZVals1 = length(ZVals1)
nZVals2 = length(ZVals2)

initcols = c('blue', 'orange')
cols = colorRampPalette(c('blue', 'orange'))(length(ZVals1)*length(ZVals2))

## Parameters fixed over entire figure
dvals <- unique(parmat$d)
gamvvals <- unique(parmat$gamv)
##gampvals <- unique(parmat$gamp)
Rpvals <- unique(parmat$Rp)
tbvals <- unique(parmat$tb)
rhovals <- unique(parmat$rho)
pmuvals <- unique(parmat$pmu)

gridvals <- expand.grid(d = dvals, tb = tbvals,
                        gamv = gamvvals, rho = rhovals,
                        pmu = pmuvals)

for(ii in 1:nrow(gridvals)){

    dval <- gridvals[ii,'d']
    tbval <- gridvals[ii,'tb']
    gamvval <- gridvals[ii,'gamv']
    gampval <- gridvals[ii,'gamp']
    rhoval <- gridvals[ii,'rho']
    pmuval <- gridvals[ii,'pmu']

    col1 = 'blue'
    col2 = 'orange'

    ## Graph both file types
    for(type in c('.eps')){

        wi.fig <- parmat$d == dval &
            parmat$gamv == gamvval &
            parmat$rho== rhoval &
            parmat$pmu == pmuval ## parameters to keep in the fig

        dn <- round(1/dval/365,1)
        gvn <- round(1/gamvval)

        dir.create(paste('Figures_1/'), showWarnings = FALSE)
        dir.create(paste('Figures_1/LS_',dn, '_rho_' , rhoval , sep = ''), showWarnings = FALSE)

        fn <- paste( 'Figures_1/', paste('LS_',dn, '_rho_' , rhoval , sep = ''), '/',
                    'gamv_', gvn, '_pmu_', pmuval,
                    '.eps', sep = '' )

        setEPS()
        postscript(file = fn, height = 3, width = 6)
        ##png(file = fn, height = 3, width = 6, res = 400, units = 'in')
        layout(matrix(c(1,2), nrow = 1, byrow = T), widths = c(0.4,0.4),
               heights = c(0.25,0.25))
        par(oma = c(0,2,0,0.5), mar = c(3, 1.5, 1.5, 0))


        for(i.f in 1:nFVals){
            FVal = FVals[i.f]
            plot.mat = matrix(nrow = nXVals, ncol = nZVals1*nZVals2)
            for(i.z1 in 1:nZVals1){
                for(i.z2 in 1:nZVals2){
                    wi <- parmat[,ZVal1Name] == ZVals1[i.z1] &
                        parmat[,ZVal2Name] == ZVals2[i.z2] &
                        parmat[,FValName] == FVals[i.f] & wi.fig
                    plot.mat[,2*(i.z1-1) + i.z2] = parmat[wi,YValName]
                }
            }


            wi.maxs <- apply(plot.mat, 2, function(x){XVals[which.max(x)]})
            max.vals <-apply(plot.mat, 2, max)
            matplot(NA, type = 'l', xlim = c(0,365), ylim = c(0.0,1.25),
                    xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
            polygon(x = c(0.0, tbval, tbval, 0.0), y = c(-2,-2,1.2,1.2), col = 'gray90',
                    lty = 3, lwd = 0.5, border = NA)
            polygon(x = c(0.0, tbval, tbval, 0.0), y = c(0.9,0.9,1.2,1.2), col = 'white',
                    lty = 3, lwd = 0.5, border = NA)
            matlines(XVals, plot.mat, lwd = 1.5, pch = c(1,1,2,2), lty = c(1,2,1,2),
                     type = 'l', col = c(col1, col1, col2, col2))
            points(wi.maxs, max.vals, pch = 19, cex = 0.5, col = 'black')
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

                legendtext = as.expression(sapply(1:length(ZVals1),
                                                  function(i){bquote('R'['0,p']*' = '*.(ZVals1[i]))}))
                legend(x = -10, y = 1.25, legend = legendtext, col = c(col1,col2),
                       lwd = 2, bty = 'n', cex = 0.7,
                       title = as.expression(bquote(underline('Pathogen Transmission'))))


                legendtext = as.expression(sapply(1:length(ZVals2),
                                                  function(i){bquote(gamma[p]*' = '*.(round(1/ZVals2[i])))}))
                legend(x = 0.55*365, y = 1.25, legend = legendtext, col = 'black',
                       lwd = 2, bty = 'n', cex = 0.7, lty = c(1,2),
                       title = as.expression(bquote(underline('Pathogen Recovery'))))

            }

        }
        dev.off()
    } ## End loop through figure types
} ## Loop through gridvals (paramters that are fixed in the figure)
