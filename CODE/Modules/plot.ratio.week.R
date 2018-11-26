####################################################################################################
# 'plot.ratio.week.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This routine plots the count ratios against the week they correspond to, for a specific organ
# found in the global environment and set by main.R. By 'week,' here, we mean the week described
# (per dataset) in weeks.csv (in the dataset folder). The count ratio per dataset is averaged over
# duplicates in that dataset (see quant.prof.R for details). If multiple datasets correspond to one
# and the same week, we additionally average those datasets as well (but still display the original,
# duplicate-averaged datasets as well).
#
# (ratios of masses held by left and right regions) for datasets of current organ [(fils.v)x(H2O,D2O) df]
#
# This routine operates as plot.ratio.R but for the data produced by quant.fix.R, instead of by
# quant.prof.R.
#
# * Plots of all datasets for the current organ, together with (H2O/D2O) pure profiles
#   and vertical lines at the (sliding; per dataset) left/right cutoffs and valley.
#
#   DEPENDENCIES
# <- settings
# <- presets
# <- all modules before it in main.R
# -> all modules after it in main.R
#
####################################################################################################
# AGGREGATE RATIOS AND AVERAGE (ABSOLUTE) SKEWNESS RATIOS PER WEEK -->>

 # Filenames of original (excl. controls) and repaired (excl. controls and pure pollutants) datasets [(orig,fixd)-list of fils.v subsets]
  fils.trim.l <- list(
    orig = dat.fils.df$fil[-which(dat.fils.df$fil %in% ctrl.fils.df$fil)] , # dataset filenames excluding controls
    fixd = dat.fils.df$fil[-which(dat.fils.df$fil %in% ctrl.fils.df$fil)] ) # filenames excluding controls and pure pollutants

 # List of lists of (L,R,T) ratios per dataset for original and repaired data [(orig,fixd)-list of lists]
  ratio.trim.l <- list(
    orig = lapply(cnt.l[fils.trim.l$orig] , # list original ratios [fils.trim.l$orig-named list of (L,R,T)x(A,B,AB) df's]
                  function(arg){out <- arg[,c('A','B','AB')] ; names(out) <- rownames(arg) ; colnames(out) <- c('A','B','AB') ; return(out)}) ,
    fixd = lapply(cnt.fix.l[fils.trim.l$fixd] , # list repaired ratios [fils.trim.l$fixd-named list of (L,R,T)x(A,B,AB) df's]
                  function(arg){out <- arg[,c('A','B','AB')] ; names(out) <- rownames(arg) ; colnames(out) <- c('A','B','AB') ; return(out)}) )

 # Update ratio.SR.week.aggr.l with (weekly, average, L/R/T) ratios and SRs [(orig,fix)-list of (week,L,R,T,SR)-col df]
  for(c.fix.stat in c('orig','fixd'))
  { # treat original/repaired data independently (iterate over repair-status)
   for(c.week in as.character(weeks.v))
   { # iterate over week numbers represented in current organ

    # Dataset filenames corresponding to current week (under current repair-status) [fils.v subset]
     c.fils <- weeks.l[[c.week]][weeks.l[[c.week]] %in% fils.trim.l[[c.fix.stat]]]
    # Ratios corresponding to current week (under current repair-status) [c.fils-list of (L,R,T)x(A,B,AB) df's]
     c.ratio.trim.l <- ratio.trim.l[[c.fix.stat]][c.fils]
    # Compile df holding mean (i.e. AB-)ratios of datasets corresponding to current week (under current repair-status) [(L,R,T)x(c.fils) df]
     aux.AB.df <- data.frame( matrix(NA,nrow=length(ratio.rows),ncol=length(c.fils)) ) # NA-preallocate df
       dimnames(aux.AB.df) <- list(ratio.rows , c.fils) # name rows after intervals and columns after dataset fillenames
     for(c.fil in c.fils){aux.AB.df[,c.fil] <- c.ratio.trim.l[[c.fil]][,'AB']} # copy AB cols (one by one) into corresponding df cols
    # Average (L,R,T) AB-ratios and (absolute, max-across-H2O/D2O-incl-duplicates) SRs over current week datasets [(week,L,R,T,SR)-col df]
     ratio.SR.week.aggr.l[[c.fix.stat]][ratio.SR.week.aggr.l[[c.fix.stat]]$week==c.week,names(rowMeans(aux.AB.df))] <- rowMeans(aux.AB.df) # ratios
     ratio.SR.week.aggr.l[[c.fix.stat]][ratio.SR.week.aggr.l[[c.fix.stat]]$week==c.week,'prp'] <- mean(unlist( abs(SR.df[c.fils,dupls]) )) # SRs

   } # end of iteration over week numbers represented in current organ
   # Remove NaN entries (i.e. weeks w/o datasets; may occur due to removal of datasets designated 'pure pollutants')
    ratio.SR.week.aggr.l[[c.fix.stat]] <- ratio.SR.week.aggr.l[[c.fix.stat]][!apply(sapply(ratio.SR.week.aggr.l[[c.fix.stat]],is.nan) , 1 , any),]
  } # end of treating original/repaired data independently

 # List holding individual [(ratio,SR)-list of (week,A,B,AB)-col df]
  for(c.fil in dat.fils.df$fil)
  { # retain ratios (per duplicate and mean) [single-row (week,A,B,AB)-col df]
   # Ratios
    aux.df <- data.frame( row.names=names(cnt.fix.l[c.fil]) ,
      week = weeks.df[weeks.df$Filename==c.fil,'Weeks'] , # dataset week
      cnt.fix.l[[c.fil]]['R',setdiff(ratio.cols,dat.cols)] ) # dataset (A,B,AB) ratios
    ratio.prp.l$ratio <- rbind(ratio.prp.l$ratio , aux.df) # attach extracted dataset ratios to overall df
   # Skewness ratios
    ratio.prp.l$prp[c.fil,'week'] <- weeks.df[weeks.df$Filename==c.fil,'Weeks'] # dataset week
    ratio.prp.l$prp[c.fil,'A'] <- max(abs(SR.df[c.fil,'A'])) # dataset duplicate A
    ratio.prp.l$prp[c.fil,'B'] <- max(abs(SR.df[c.fil,'B'])) # dataset duplicate B
    ratio.prp.l$prp[c.fil,'AB'] <- rowMeans(ratio.prp.l$prp[c.fil,dupls]) # mean across dataset duplicates
  }
  # Trim list elements by discarding controls and pure pollutants
   ratio.prp.l$ratio <- ratio.prp.l$ratio[fils.trim.l[['fixd']] , ]
   ratio.prp.l$prp <- ratio.prp.l$prp[fils.trim.l[['fixd']] , ]


# AGGREGATE RATIOS AND AVERAGE (ABSOLUTE) SKEWNESS RATIOS PER WEEK <<--
####################################################################################################
# PLOT RATIOS VS WEEKS -->>

   pdf(paste('../Results/',organ,'/all_ratios_vs_weeks_',organ,'.pdf',sep='')) # open plotting connection

 # Plot settings
  # SR.max <- -Inf # initialize max D2O SR-value across datasets
  # for(c.organ in organs) # compute max D2O SR-value across datasets [positive scalar]
  #   {SR.max <- max(SR.max , max(abs( SRs.l[[c.organ]][!(rownames(SRs.l[[c.organ]]) %in% names(fils.poll.l[[c.organ]])) , cols$D2O] )) )}
  apply(sapply(ratio.SR.week.aggr.l[['fixd']] , is.nan) , 1 , any)
  plot.max.y <- (1+plot.eps.df$y)*max(sapply(ratio.SR.week.aggr.l,function(arg){max(arg$R)})) # max y-range [scalar]
    plot.N.col <- 1E3 # no.of colors in palette
  plot.col.pal <- colorRampPalette(c('black','white'))(plot.N.col) # make grayscale palette [plot.N.col-long string vector]

 # Plot mean ratios
  aux.rats <- ratio.SR.week.aggr.l[['orig']] # aggregated data (original)
    plot.col.dat <- plot.col.pal[as.integer(1+(plot.N.col-1)*round(aux.rats$prp/max(aux.rats$prp),3))] # colors per prp-value [plot.N.col-long nat vector]
  aux.rats <- ratio.SR.week.aggr.l[['fixd']] # aggregated data (repaired)
    plot.col.dat <- plot.col.pal[as.integer(1+(plot.N.col-1)*(1-round(aux.rats$prp/max(aux.rats$prp),3)))] # colors per prp-value [plot.N.col-long nat vector]
  plot(aux.rats$week , aux.rats$R , typ='p' , pch=1 , cex=2.4 , col='red' , 
    ylim=c(0,0.5) , xlab='week' , ylab='deuteration' , main=organ) # hollow points (repaired)
  lines(aux.rats$week , aux.rats$R , typ='p' , pch=19 , cex=2 , col=plot.col.dat) # points (repaired)
  lines(aux.rats$week , aux.rats$R , typ='l' , pch=16 , col='red') # joining lines (repaired)

 # Plot individual ratios (incl duplicates)
  for(c.week in as.character(weeks.v))
  {
   # Dataset filenames corresponding to current week that have been retained (under current repair-status) [fils.v subset]
    c.fils <- weeks.l[[c.week]][weeks.l[[c.week]] %in% fils.trim.l[[c.fix.stat]]]
   # Ratios corresponding to current week (under current repair-status) [c.fils-list of (L,R,T)x(A,B,AB) df's]
    c.ratio.trim.l <- ratio.trim.l[[c.fix.stat]][c.fils]
   # Compile df holding AB-ratios of all datasets corresponding to current week (under current repair-status) [(L,R,T)xc.fils df]
    aux.AB.df <- data.frame( matrix(NA,nrow=length(ratio.rows),ncol=length(c.fils)) ) # NA-preallocate df
      dimnames(aux.AB.df) <- list(ratio.rows , c.fils) # name rows after intervals and columns after dataset fillenames
    for(c.fil in c.fils){aux.AB.df[,c.fil] <- c.ratio.trim.l[[c.fil]][,'AB']} # copy AB cols (one by one) into corresponding df cols
   # Plot individual datasets
    lines(rep(as.integer(c.week),dim(aux.AB.df)[2]) , aux.AB.df['R',] , typ='p' , col='blue' , pch=8 , cex=0.5)
  }

 # Plot individual ratios (incl duplicates)
  plot.col.dat.l <- list(
    A = plot.col.pal[as.integer(1+(plot.N.col-1)*(1-round(ratio.prp.l$prp$A/max(ratio.prp.l$prp[,dupls]),3)))] ,
    B = plot.col.pal[as.integer(1+(plot.N.col-1)*(1-round(ratio.prp.l$prp$B/max(ratio.prp.l$prp[,dupls]),3)))] )
  names(plot.col.dat.l$A) <- rownames(ratio.prp.l$prp)
  names(plot.col.dat.l$B) <- rownames(ratio.prp.l$prp)
  for(c.fil in fils.trim.l[['fixd']])
  { # iterate over datasets (excluding contorls and pure pollutants)
   # Plot individual datasets
    jttr <- (runif(1)-0.5)*data.frame(A=1,B=-1) # jitter (for plotting points so as not to overlap)
    lines(jttr+rep(ratio.prp.l$ratio[c.fil,]$week,2) , ratio.prp.l$ratio[c.fil,dupls] , typ='l' , col='black')
    lines(jttr$A+ratio.prp.l$ratio[c.fil,]$week , ratio.prp.l$ratio[c.fil,'A'] , typ='p' , col=plot.col.dat.l$A[[c.fil]] , pch=19 , cex=1)
    lines(jttr$B+ratio.prp.l$ratio[c.fil,]$week , ratio.prp.l$ratio[c.fil,'B'] , typ='p' , col=plot.col.dat.l$B[[c.fil]] , pch=19 , cex=1)
  }

   dev.off() # close plotting connection

# PLOT RATIOS VS WEEKS <<--
####################################################################################################