####################################################################################################
# 'pure.locs.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This routine identifies peaks, valleys and cutoffs for the pure profiles of the organ set in
# settings.R. Here also, we distinguish start/end, H2O/D2O pure profiles; cf. pure.tab.R.
#
# For each position (start/end) and molecular type (H2O/D2O), we start by assigning the pure profile
# from the corresponding element of pure.prof.l to c.dat.df. The two peaks and the valley are then
# found using the routine extrema.locs.R, see its preamble for details. Cutoffs, on the other hand,
# are timepoints bounding 100*(1-peak.thr)% of the distribution mass (i.e. AUC) of the pure profile
# (that AUC is actually one, cf. pure.tab.R). The cutoffs depend on how we divide the remaining
# 100*peak.thr% of the AUC (to-be-discarded) to left/right tails. Here, we divide it proportionally
# to the AUC held by the (uncropped) left and right half-peaks.
#
# Specifically, we store the AUCs held by left, middle and right intervals (i.e. left of the left
# peak, in between the peaks and right of the right peak) in the (L,M,R)-vector aux.AUC; their sum
# is one. The bounding interval certainly contains the M (middle) AUC (equaling aux.AUC$M), so it
# remains to split what's left of the AUC to the left/right (L/R) halves of the L/R peaks. The AUC
# to be split equals 1 - peak.thr - aux.AUC$M, and it is it split proportionally (see above), i.e.
# according to the ratio aux.AUC$L/aux.AUC$R. The cutoffs are found by comparing the cumulative sums
# for the L/R regions to the targets (stored in the (L,R)-named vector aux.AUC.trg). Results are
# stored in the corresponding entries of pure.loc.l.
#
#   INPUTS
# None.
#
#   OUTPUTS
# * pure.loc.l = two-element list of lists, w/ each element a list having an H2O and a D2O elements;
#                these elements are (x,idx)-columned df's recording the left/right peak locations
#                (rows Lmax anr Rmax) and left/right boundaries for the (averaged, mollified, H2O
#                and D2O) controls of the current position (start/end) and current organ
#                            [(start,end)-list of (H2O,D2O)-lists of (L,Lmax,M,Rmax,R)x(idx,x) df's]
#
#   DEPENDENCIES
# <- settings.R
# <- presets.R
# <- pure.tab.R
# <- extrema.locs.R
# -> all modules after it in main.R
####################################################################################################
for(c.pos in poss)
{ # iterate analysis over start/end pure profiles
 for(c.mol in mols)
 { # iterate analysis over undeuterated/deuterated pure profiles

  # Load pure profile corresponding to current organ/position/type and specify data col as y-col [(N.t)x(idx,x,y) df]
   c.dat.df <- setNames(pure.prof.l[[c.pos]][,c('idx','x',c.mol)] , c('idx','x','y'))

  # Extract and record left/right peaks, as well as valley falling between them
   aux.l <- extrema.locs(c.dat.df,smoothen=FALSE) # peaks and inter-peak valley [(max,min)-list of (idx,x)-df's]
     if(dim(aux.l$max)[1]!=2){stop('ABORT: more/less than two peaks found [pure.prof.locs]')} # sanity check
   pure.loc.l[[c.pos]][[c.mol]][c('Lmax','Rmax'),] <- aux.l$max # record peak locations in pure.loc.l
   aux.l$min <- subset(aux.l$min , min(aux.l$max$x) < x & x < max(aux.l$max$x)) # retain inter-peak valleys only
     if(dim(aux.l$min)[1]!=1){stop('ABORT: no single inter-peak valley found [pure.prof.locs]')} # sanity check
   pure.loc.l[[c.pos]][[c.mol]]['M',] <- aux.l$min # record valley location in pure.loc.l

  # Extract and record boundaries
   aux.AUC <- dt.g*data.frame( # AUCs below left peak, between the two peaks, and above right peak [(L,M,R) num vector]
     L = sum(c.dat.df[c.dat.df$x < min(aux.l$max$x) , 'y']) ,
     M = sum(c.dat.df[min(aux.l$max$x) <= c.dat.df$x & c.dat.df$x <= max(aux.l$max$x) , 'y']) ,
     R = sum(c.dat.df[max(aux.l$max$x) < c.dat.df$x , 'y']) )
   if(peak.thr<=1-aux.AUC$M)
   { # left/right target AUCs determining cutoffs [(L,R) num vector]
    aux.AUC.trg <- (1-peak.thr-aux.AUC$M)/(1-aux.AUC$M)*data.frame('L'=aux.AUC$L , 'R'=aux.AUC$R)
   }else{ stop(sprintf('ABORT: AUC to discard (%.02f) exceeds tail AUC (%.02f) [pure.prof.locs]',peak.thr,1-aux.AUC$M)) }
   # Find and record in pure.loc.l rightmost/leftmost gridpoints hitting left/right target
    c.crop.l <- list(
      R = c.dat.df[max(aux.l$max$x) < c.dat.df$x,] , # df of entries with times above right peak [sub-df of c.dat.df]
      L = data.frame(lapply(c.dat.df[c.dat.df$x < min(aux.l$max$x),],rev)) ) # df of entries with times below left peak [upside-down sub-df of c.dat.df]
    aux.idx <- list(
      R = max(c.crop.l$R[dt.g*cumsum(c.crop.l$R$y) < aux.AUC.trg$R , 'idx']) , # index of right cutoff [entry of c.dat.df$idx]
      L = min(c.crop.l$L[dt.g*cumsum(c.crop.l$L$y) < aux.AUC.trg$L , 'idx']) ) # index of left cutoff [entry of c.dat.df$idx]
    pure.loc.l[[c.pos]][[c.mol]][c('L','R'),] <- # record left/right cutoffs
      rbind(c.dat.df[c.dat.df$idx==aux.idx$L , c('idx','x')] , c.dat.df[c.dat.df$idx==aux.idx$R , c('idx','x')])

  # Check that ordering is correct
   if(any(pure.loc.l[[c.pos]][[c.mol]][c('L','Lmax','M','Rmax','R'),'x'] != cummax(pure.loc.l[[c.pos]][[c.mol]][c('L','Lmax','M','Rmax','R'),'x'])))
    {stop('ABORT: bounds-peaks-valley incorrectly ordered [pure.locs.R]')}

 } # end of iterating analysis over undeuterated/deuterated pure profiles
} # end of iterating analysis over start/end pure profiles

    print(sprintf('Feature identification of pure profiles complete [%s]',organ))