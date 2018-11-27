####################################################################################################
# 'pure.tab.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This routine tabulates pure start/end profiles for the organ (specified in settings.R) over the
# uniform timegrid t.g (see presets.R). By virtue of the experimental design, controls are measured
# (typically) at the start and end of the GC/MS run in which all samples from a specific organ are
# processed. This is evident in the weeks.csv file (in the dataset folder), as control files (marked
# as 'CTRL') occupy the top/bottom rows. (In the few exceptions where either start or end controls
# are missing, we simply duplicate the end/start controls for our analysis here.)
#
# We treat start/end control ensembles independently. (Below, we refer to the 'start/end' label as
# 'position.') For each position (i.e. for each such ensemble), we form the list c.ctrl.df.l with
# elements H2O and D2O (referred throughout as 'molecular types'). Each such element is a df with an
# index column ('idx'), a time grid column ('x') and twice (because of duplicate measurements) as
# many extra columns as control datasets at that position. For given position and molecular type, we
# record in these extra columns the control counts of that molecular type and position. Each column
# is then normalized so that its area-under-the-curve (AUC=sum of column entries multiplied by the
# timestep dt.g - see presets.R) equals one; this makes each column a tabulated distribution. The
# average of these distributions is also a distribution, corresponding to the average control
# profile for that given position and molecular type. In this manner, we obtain four distros, which
# we then mollify (using a Gaussian kernel, see presets.R) and renormalize to their original AUC
# (i.e. to one). These final distributions are saved as elements of the list pure.prof.l, which is a
# start/end-list whose elements (df's) have columns 'idx', 'x', 'H2O' and 'D2O'; these obviously
# tabulate the pure H2O/D2O profiles derived from the start/end control ensembles.
#
#   INPUTS
# None.
#
#   OUTPUTS
# * pure.prof.l = (H2O,D2O)-list each element of which is a df holding an index column, a timegrid
#                 column, an H2O column tabulating the (average, mollified) undeuterated pure
#                 profile over the timegrid and a D2O column tabulating the (average, mollified)
#                 deuterated pure profile over the same timegrid
#                                                   [(start,end)-list of (N.t)x(idx,x,H2O,D2O) df's]
#
#   DEPENDENCIES
# <- settings.R
# <- presets.R
# -> all modules after it in main.R
####################################################################################################
for(c.pos in poss)
{ # treat start/end profiles separately

 # Initialize (H2O,D2O)-list to hold data for current position (extra columns added below) [(H2O,D2O)-list of (idx,x)-col df's]
  c.ctrl.df.l <- list(
    H2O = data.frame(idx=idx.g , x=t.g) , # initialize undeuterated control data [df]
    D2O = data.frame(idx=idx.g , x=t.g) ) # initialize deuterated control data [df]

 # List df's w/ columns holding index, timegrid and data (normalized to unit AUC) for current position [(H2O,D2O)-list of df's]
  for(c.fil in ctrl.fils.split.l[[c.pos]]$fil)
  { # iterate over control datasets belonging to current position (stored in corresponding element of ctrl.fils.split.l)
   # Load the current control dataset (truncated and repaired, cf. presets.R) [(N.t)x(dat.cols) df]
    c.dat.df <- dat.tab.l[[c.fil]]
   # Total AUCs for current control dataset (in duplicate) per molecular type [(H20,D20)-named list of two-element num vectors]
    aux.dat.norms <- list(H2O=dt.g*colSums(c.dat.df[,cols$H2O]) , D2O=dt.g*colSums(c.dat.df[,cols$D2O]))
   # List normalized profiles (concatenated with index vector and time grid) [(H20,D20)-named list of df's]
    c.ctrl.df.l <- list(
      H2O = cbind(c.ctrl.df.l$H2O , sweep(c.dat.df[,cols$H2O],2,aux.dat.norms$H2O,'/')) , # normalize and attach current H2O control data [df]
      D2O = cbind(c.ctrl.df.l$D2O , sweep(c.dat.df[,cols$D2O],2,aux.dat.norms$D2O,'/')) ) # normalize and attach current D2O control data [df]
  } # end of iteration over control datasets

 # Extract and list in pure.prof.l the pure (start/end, H2O/D2O) profiles [(start,end)-list of (idx,x,H2O,D2O) df's]
  aux.avg.cols <-  # indices of data-holding cols in c.ctrl.df.l element [nat vector]
  pure.prof.l[[c.pos]] <- data.frame( # average H2O/D2O data (incl. duplicates) and attach to index and time grid
    idx = idx.g , # as per the initialization of the list above
    x = t.g , # as per the initialization of the list above
    H2O = rowMeans(c.ctrl.df.l$H2O[ , which(!(colnames(c.ctrl.df.l$H2O) %in% c('idx','x'))) ]) , # across-data-column means
    D2O = rowMeans(c.ctrl.df.l$D2O[ , which(!(colnames(c.ctrl.df.l$D2O) %in% c('idx','x'))) ]) ) # across-data-column means
  for(c.mol in mols)
  { # treat (mollify, renormalize and remove offset from) undeuterated/deuterated profiles
   pure.prof.l[[c.pos]][,c.mol] <- mollify( data.frame('x'=pure.prof.l[[c.pos]]$x , 'y'=pure.prof.l[[c.pos]][,c.mol]) )$y
   pure.prof.l[[c.pos]][,c.mol] <- pure.prof.l[[c.pos]][,c.mol] - min(pure.prof.l[[c.pos]][,c.mol]) # remove offset
  } # end of treating undeuterated/deuterated profiles

} # end separate treatment of start/end profiles

    print(sprintf('Tabulation of pure profiles complete [%s]',organ))