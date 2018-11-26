####################################################################################################
# 'plot.profs.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This routine plots the count profiles for all datasets corresponding to the organ current in the
# global environment (set by main.R). The profiles are normalized to have unit AUC over the rough
# interval (cf. settings.R) and are only plotted  over it. An important part of this module is
# setting the "smart cutoff intervals" per dataset, i.e. deciding where to set the lower and upper bounds
# (as well as the valley) for each dataset in the current organ. The way this is done is by linear
# interpolation between (the corresponding times for) start and end profiles, with the datasets
# ordered as in weeks.df. These "sliding" cutoffs are also included (as vertical lines) in the plot,
# but all portrayed profiles (datasets, pure profiles, pollutant profiles) are normalized over the
# rough interval and not over the corresponding "sliding" intervals.
#
#   INPUTS
# None.
#
#   OUTPUTS
# * Plots of all (H2O/D2O) datasets (blue/red) for the current organ, together with (H2O/D2O) pure
#   profiles (light blue/pink), the pollutant profile (green) and vertical lines at the (sliding;  
#   per dataset) left/right (H2O/D2O) cutoffs and valley (blue/red).
#
#   DEPENDENCIES
# <- settings.R
# <- presets.R
# -> main.R
####################################################################################################
# Plot settings
 plot.eps.df <- data.frame(x=2E-3,y=2E-2) # plotting buffers [single row, (x,y)-col df]
 plot.max.y <- (1+plot.eps.df$y)*max(pure.prof.l$start[,c('H2O','D2O')] , pure.prof.l$end[,c('H2O','D2O')]) # max y [scalar]
 plot.label.lims <- list( # start/end coordinates for inserting 'start', 'end' plot labels [(start,end)-list of single-row (x,y)-df's]
   start = data.frame( x=t.g[which.max(pure.prof.l$start$H2O)] , y=(1+plot.eps.df$y)*max(pure.prof.l$start$H2O) ) ,
   end = data.frame( x=t.g[which.max(pure.prof.l$end$H2O)] , y=(1+plot.eps.df$y)*max(pure.prof.l$end$H2O) ) )
 plot.label.lims$start$x <- (1 + 2*((plot.label.lims$start$x>plot.label.lims$end$x) - 1/2)*plot.eps.df$x)*plot.label.lims$start$x # 'start' label position
 plot.label.lims$end$x <- (1 - 2*((plot.label.lims$start$x>plot.label.lims$end$x) - 1/2)*plot.eps.df$x)*plot.label.lims$end$x # 'end' label position

# SR-plots for all non-pollutant datasets of current organ
 if(plot.swtc)
 { # check plot switch
    pdf(paste('../Results/',organ,'/SR_plots_',organ,'.pdf',sep='')) # open plotting connection
  c.df <- SR.df[sort(rowMeans(SR.df[,cols$D2O]) , index.return=T)$ix , ] # order df by increasing SR of D2O datasets (averaged over duplicates)
  layout(matrix(c(1,2), nrow=1 , ncol=2)) # arrange plots in 1x2 matrix
  # Plot H2O (blue) and D2O (red) SRs of duplicates (x-axis: A-duplicate; y-axis: B-duplicate)
   plot(c.df[,cols$D2O] , typ='n' , # initialize plot (with right bounds) without plotting anything
     xlim=range(c.df[,cols$D2O[1]],c.df[,cols$H2O[1]]) , ylim=range(c.df[,cols$D2O[2]],c.df[,cols$H2O[2]]) ,
     xlab='A' , ylab='B' , main='A and B SRs (H2O, D2O)' )
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="SlateGray4") # background color
   lines(c.df[rownames(c.df) %in% ctrl.fils.df$fil,][,cols$D2O] , typ='p' , col=alpha('HotPink1',0.25) , pch=16 , cex=1) # D2O controls
   lines(c.df[!rownames(c.df) %in% ctrl.fils.df$fil,][,cols$D2O] , typ='p' , col='HotPink1' , pch=16 , cex=1) # D2O datasets
   lines(c.df[rownames(c.df) %in% ctrl.fils.df$fil,][,cols$H2O] , typ='p' , col=alpha('LightSkyBlue1',0.25) , pch=16 , cex=0.5) # H2O controls
   lines(c.df[!rownames(c.df) %in% ctrl.fils.df$fil,][,cols$H2O] , typ='p' , col='LightSkyBlue1' , pch=16 , cex=0.5) # H2O datasets
  # Plot A (green) and B (yellow) SRs of H2O-and-D2O datasets (x-axis: H2O-values; y-axis: D2O-values)
   plot(c.df[,cols$D2O] , typ='n' ,  # initialize plot (with right bounds) without plotting anything
     xlim=range(c.df[,cols$A[1]],c.df[,cols$B[1]]) , ylim=range(c.df[,cols$A[2]],c.df[,cols$B[2]]) ,
     xlab='H2O' , ylab='D2O' , main='H2O and D2O SRs (A, B)' )
     rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="SlateGray4") # background color
   lines(c.df[rownames(c.df) %in% ctrl.fils.df$fil,][,cols$A] , typ='p' , col=alpha('LimeGreen',0.25) , pch=16 , cex=1) # A-duplicate controls
   lines(c.df[!rownames(c.df) %in% ctrl.fils.df$fil,][,cols$A] , typ='p' , col='LimeGreen' , pch=16 , cex=1) # A-duplicate datasets
   lines(c.df[rownames(c.df) %in% ctrl.fils.df$fil,][,cols$B] , typ='p' , col=alpha('Tan1',0.25) , pch=16 , cex=1) # B-duplicate controls
   lines(c.df[!rownames(c.df) %in% ctrl.fils.df$fil,][,cols$B] , typ='p' , col='Tan1' , pch=16 , cex=1) # B-duplicate datasets
   for(c.row in rownames(c.df)) # join duplicates (A- and B-datasets)
   {
    if(c.row %in% ctrl.fils.df$fil)
      {lines(c(c.df[c.row,cols$A[1]],c.df[c.row,cols$B[1]]) , c(c.df[c.row,cols$A[2]],c.df[c.row,cols$B[2]]) , typ='l' , col=alpha('MintCream',0.25))}
    if(!c.row %in% ctrl.fils.df$fil)
      {lines(c(c.df[c.row,cols$A[1]],c.df[c.row,cols$B[1]]) , c(c.df[c.row,cols$A[2]],c.df[c.row,cols$B[2]]) , typ='l' , col='MintCream')}
   }
     dev.off() # close plotting connection
 } # end of check plot switch

# 2x2 combi-plots for all non-pollutant datasets of current organ
 for(c.fil in dat.fils.df$fil)
 {
  # If plot switch is on, plot datasets together with peaks/valleys/cutoffs/pure- and valley-derived candidate cutoffs
   if(plot.swtc)
   { # check plot switch    
       pdf(plot.fils[c.fil,'combi']) # open plotting connection
     layout(matrix(c(1,2,3,4), nrow=2 , ncol=2 , byrow = FALSE))
     for(c.dupl in dupls)
     { # iterate over dataset duplicates
      # Load profiles for current dataset [(N.t)x(idx,x,dat.cols) df]
       c.cnt.df <- dat.tab.l[[c.fil]]
      # Store (H2O/D2O) count-profiles for current duplicate [(N.t)x(idx,x,H2O,D2O) df]
       aux.df <- c.cnt.df[ , c('idx','x',cols[[c.dupl]])] # excise non-current duplicate
         for(c.col in cols[[c.dupl]]){names(aux.df)[names(aux.df)==c.col] <- strsplit(c.col,'[.]')[[1]][1]} # rename data cols
      # Store extrema of (H2O/D2O) count-profiles for current duplicate [(H2O,D2O)-list of (max,min)-lists of (loc,idx)-lists]
       c.extr.l <- lapply(mols , function(arg){extrema.locs(cbind(aux.df[,c('idx','x')],y=aux.df[,arg]),smoothen=TRUE)})
         names(c.extr.l) <- mols
      # Comparison between current duplicate and start/end pure profiles [mols-list of (start,end)x(idx,x,shft,mad) df's]
       c.peak.l <- dat.vs.pure(aux.df)
      for(c.mol in mols)
      { # iterate over molecular types (H2O/D2O) of current duplicate
       # Element of dat.cols corresponding to current molecular type and duplicate [string]
        c.dat.col <- intersect(cols[[c.mol]],cols[[c.dupl]])
       # Mollified version of part of aux.df corresponding to current molecular type [(N.t)x(x,y) df]
        aux.moll.df <- mollify(cbind(aux.df[,c('idx','x')],y=aux.df[,c.mol]))[,c('x','y')]
       # Plots
        plot(cbind(x=aux.df$x , y=aux.df[,c.mol]) , typ='l' , col='grey' , main=c.dat.col) # prof
        lines(aux.moll.df$x , aux.moll.df$y , typ='l') # mollified prof
        lines(rep(prof.loc.l[[c.fil]]$x['R',c.dat.col],2) , c(0,par("usr")[4]) , typ='l' , col='Yellow1' , lwd=3) # R-cut
        lines(rep(prof.loc.l[[c.fil]]$x['L',c.dat.col],2) , c(0,par("usr")[4]) , typ='l' , col='Yellow1' , lwd=3) # L-cut
        lines(rep(prof.pure.cand.l[[c.fil]]$x['R',c.dat.col] , 2) , c(0 , par("usr")[4]) , typ='l' , col='OliveDrab3') # pure R-line
        lines(rep(prof.pure.cand.l[[c.fil]]$x['L',c.dat.col] , 2) , c(0 , par("usr")[4]) , typ='l' , col='OliveDrab1') # pure L-line
        lines(rep(prof.valley.cand.l[[c.fil]]$x['R',c.dat.col] , 2) , c(0 , par("usr")[4]) , typ='l' , col='Plum3') # valley R-line
        lines(rep(prof.valley.cand.l[[c.fil]]$x['L',c.dat.col] , 2) , c(0 , par("usr")[4]) , typ='l' , col='Plum1') # valley L-line
        lines(c.extr.l[[c.mol]]$max$x , aux.df[aux.df$idx %in% c.extr.l[[c.mol]]$max$idx,c.mol] , typ='p' , pch=16 , col='blue') # peaks
        lines(c.extr.l[[c.mol]]$min$x , aux.df[aux.df$idx %in% c.extr.l[[c.mol]]$min$idx,c.mol] , typ='p' , pch=16 , col='red') # valleys
        lines(rep(c.peak.l[[c.mol]]$pp.x[1],2) , c(0,aux.df[aux.df$idx %in% c.peak.l[[c.mol]]$pp.idx[1],c.mol]) , typ='l' , col='blue') # PP
        lines(rep(c.peak.l[[c.mol]]$sp.x[1],2) , c(0,aux.df[aux.df$idx %in% c.peak.l[[c.mol]]$sp.idx[1],c.mol]) , typ='l' , col='lightblue') # SP
        lines(rep(prof.loc.l[[c.fil]]$x['M',c.dat.col],2) , c(0 , par("usr")[4]) , typ='l' , col='red') # midpoint line
      } # iterate over molecular types (H2O/D2O) of current duplicate
     } # end of iteration over dataset duplicates
       dev.off() # close plotting connection
   } # end of check plot switch
 }

# Profile plots (against a backdrop of start/end pure profiles and pollutant profile)
 for(c.fil in fils.v)
 {
  # Initialize (H2O,D2O)-list that will hold current dataset [(H2O,D2O)-list of (idx,x)-col df's]
   c.dat.l <- list(H2O=data.frame(idx=1:N.t , x=t.g) , D2O=data.frame(idx=1:N.t , x=t.g))
  # Normalized H2O/D2O profiles paired with index vector and timegrid [(H20,D20)-list of (N.t)x(idx,x,H2O.A,H2O.B) and (N.t)x(idx,x,D2O.A,D2O.B) df's]
   c.dat.df <- dat.tab.l[[c.fil]] # dataset corresponding to c.fil, truncated to fit the rough interval [N.t-long df]
   aux.AUC <- list(H2O=colSums(c.dat.df[,cols$H2O])*dt.g , D2O=colSums(c.dat.df[,cols$D2O])*dt.g ) # rough interval AUCs [(H2O,D2O)-list of 1x2 vectors]
   c.dat.l <- list( # prepare data for plotting [(H2O,D2O)-list of (idx,x,cols$H2O) and (idx,x,cols$D2O) df's]
     H2O = cbind(c.dat.l$H2O , sweep(c.dat.df[,cols$H2O],2,aux.AUC$H2O,'/')) , # normalize and attach current H2O control data [df]
     D2O = cbind(c.dat.l$D2O , sweep(c.dat.df[,cols$D2O],2,aux.AUC$D2O,'/')) ) # normalize and attach current D2O control data [df]
  # Load left/mid/right cutoff x-values per duplicate of current dataset [(L,M,R)x(dat.cols) df]
   c.LMR <- prof.loc.l[[c.fil]]$x
  # Plot settings
   plot.dat.cols <- list( # names of data cols [(H2O,D2O)-list of string vectors]
     H2O = colnames(c.dat.l$H2O)[which(!(colnames(c.dat.l$H2O) %in% c('idx','x')))] ,
     D2O = colnames(c.dat.l$D2O)[!(colnames(c.dat.l$D2O) %in% c('idx','x'))] )
   plot.max.y <- (1+plot.eps.df$y)*max(cbind( # max y-range [scalar]
     c.dat.l$H2O[,plot.dat.cols$H2O] , c.dat.l$D2O[,plot.dat.cols$D2O] , # data y-cols
     pure.prof.l$start[,c('H2O','D2O')] , pure.prof.l$end[,c('H2O','D2O')] )) # pure profile y-cols
  # If plot switch is on, plot H2O and D2O distributions (in duplicate) against start/end profiles
   if(plot.swtc)
   { # check plot switch
      pdf(plot.fils[c.fil,'plot']) # open plotting connection
    # Plot normalized H2O/D2O profiles (normalized over rough interval)
     plot(t.g , c.dat.l$H2O[,plot.dat.cols$H2O[1]] , typ='l' , col='blue' ,
       ylim=c(0,plot.max.y) , xlab='time' , ylab='density' , main=c.fil) # H2O profile
     lines(t.g , c.dat.l$H2O[,plot.dat.cols$H2O[2]] , typ='l' , col='blue') # H2O profile (duplicate)
     lines(t.g , c.dat.l$D2O[,plot.dat.cols$D2O[1]] , typ='l' , col='red') # D2O profile
     lines(t.g , c.dat.l$D2O[,plot.dat.cols$D2O[2]] , typ='l' , col='red') # D2O profile (duplicate)
    # Plot pure profiles (normalized over rough interval)
     lines(t.g , pure.prof.l$start$H2O , typ='l' , col='lightblue' , ylim=c(0,plot.max.y)) # start H2O profile
     lines(t.g , pure.prof.l$end$H2O , typ='l' , col='lightblue') # end H2O profile
     lines(t.g , pure.prof.l$start$D2O , typ='l' , col='pink') # start D2O profile
     lines(t.g , pure.prof.l$end$D2O , typ='l' , col='pink') # end D2O profile
     text(plot.label.lims$start$x , plot.label.lims$start$y , 'start' , col='grey') # plot label next to start profile
     text(plot.label.lims$end$x , plot.label.lims$end$y , 'end' , col='grey') # plot label next to end profile
    # Plot pollutant profile (normalized over rough interval)
     lines(poll.prof.df$x , poll.prof.df$y , typ='l' , col='lightgreen')
    # Plot vertical cutoffs
     lines(c(c.LMR['L','H2O.A'] , c.LMR['L','H2O.A']) , c(0,plot.max.y) , typ='l' , col='RoyalBlue3') # H2O left cutoff (A)
     lines(c(c.LMR['M','H2O.A'] , c.LMR['M','H2O.A']) , c(0,plot.max.y) , typ='l' , col='RoyalBlue3') # H2O mid cutoff (A)
     lines(c(c.LMR['R','H2O.A'] , c.LMR['R','H2O.A']) , c(0,plot.max.y) , typ='l' , col='RoyalBlue3') # H2O right cutoff (A)
     lines(c(c.LMR['L','H2O.B'] , c.LMR['L','H2O.B']) , c(0,plot.max.y) , typ='l' , col='RoyalBlue1') # H2O left cutoff (B)
     lines(c(c.LMR['M','H2O.B'] , c.LMR['M','H2O.B']) , c(0,plot.max.y) , typ='l' , col='RoyalBlue1') # H2O mid cutoff (B)
     lines(c(c.LMR['R','H2O.B'] , c.LMR['R','H2O.B']) , c(0,plot.max.y) , typ='l' , col='RoyalBlue1') # H2O right cutoff (B)
     lines(c(c.LMR['L','D2O.A'] , c.LMR['L','D2O.A']) , c(0,plot.max.y) , typ='l' , col='HotPink3') # D2O left cutoff (A)
     lines(c(c.LMR['M','D2O.A'] , c.LMR['M','D2O.A']) , c(0,plot.max.y) , typ='l' , col='HotPink3') # D2O mid cutoff (A)
     lines(c(c.LMR['R','D2O.B'] , c.LMR['R','D2O.B']) , c(0,plot.max.y) , typ='l' , col='HotPink3') # D2O right cutoff (A)
     lines(c(c.LMR['L','D2O.B'] , c.LMR['L','D2O.B']) , c(0,plot.max.y) , typ='l' , col='HotPink1') # D2O left cutoff (B)
     lines(c(c.LMR['M','D2O.B'] , c.LMR['M','D2O.B']) , c(0,plot.max.y) , typ='l' , col='HotPink1') # D2O mid cutoff (B)
     lines(c(c.LMR['R','D2O.B'] , c.LMR['R','D2O.B']) , c(0,plot.max.y) , typ='l' , col='HotPink1') # D2O right cutoff (B)
      dev.off() # close plotting connection
   } # end of check plot switch
 } # end of iteration over datasets for current organ