####################################################################################################
# 'dat.vs.pure.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This module accepts as input a dataset containing both H2O and D2O profiles and compares it to the
# start/end pure profiles for the organ set in settings.R and currently in the workspace.
# 
# First, the module identifies the secondary (SP) and principal (PP) peaks for the H2O and D2O
# profiles, i.e. the signal peaks corresponding to the pure profiles' left and right peaks. Contrary
# to pure profiles, dataset profile PPs need neither be rightmost (due to impurities) nor hold most
# of the AUC, due to pollutant contributing assymetrically to the secondary peak. Instead,
# we identify PPs as the peaks closest to pure profile PPs, verifying that the same peak is closest
# to both start and end pure PPs (sanity check). The SP is identified as the leftmost dataset peak,
# necessitating that the rough interval be tight enough, at its left side, to exlude spurious peaks.
#
# Once the H2O/D2O PPs are identified, each profile is shifted until its PP aligns with that of the
# start/end pure profile (considered sequentially). The pure profile is then cropped to the interval
# between its midpoint (M) and right cutoff (Rmax), has a constant (offset) subtracted bringing its
# minimum to zero and is renormalized to unit AUC. The shifted profile also undergoes the same
# procedure, and the mean average deviation between them is computed as a metric of similarity.
#
# The results are packaged in a list having one element per molecular type (H2O/D2O). That element
# is a df with the two rows corresponding to the two pure profiles and the columns detailing various
# properties of the profile: SP and PP locations, both as timegrid indices and x-values, shifts
# (from dataset profile to pure profile, measured in timegrid indices - i.e. in dt.g units) and
# finally MAD-values.
#
#   INPUTS
# * prof.in = a dataset in the form of a df, containing columns marked 'idx', 'x', 'H2O' and 'D2O'
#             and with the last two listing the H2O and D2O counts (or any renormalization thereof)
#             of the dataset [(idx,x,H2O,D2O,...)-col df]
#
#   OUTPUTS
# * out.l = list with two elements (H2O/D2O), each of which is a df holding SP index and location,
#           PP index and location, shift between (H2O/D2O) profile and (start/end) pure profiles,
#           and mean average deviation (MAD) between aligned (H2O/D2O) profile and (start/end) pure
#           profiles; the SP and PP indices and locations must be identical for both start/end
#           profiles, whereas shift and MAD will differ
#                         [(H2O/D2O)-list of (start,end)x(sp.idx,sp.x,pp.idx,pp.x,pp.shft,mad) df's]
#
#
#   DEPENDENCIES
# <- settings.R
# <- presets.R
# <- pure.profs.tab.R
# <- pure.profs.tab.R
# <- extrema.locs.R
# -> all modules after it in main.R
####################################################################################################
dat.vs.pure <- function(prof.in){

 #################################################################################################
 # Prelims
    if(!all(colnames(aux.df) %in% c('idx','x',mols))){stop('ABORT: check input profile column names [dat.vs.pure.R]')}
  # Preallocate list to hold data of input duplicate [(H2O,D2O)-list of (start,end)x(sp.idx,sp.x,pp.idx,pp.x,pp.shft,mad) df's]
   out.l <- list( # for indices of leftmost and principal peaks
     H2O = data.frame(sp.idx=rep(NA,2),sp.x=rep(NA,2),pp.idx=rep(NA,2),pp.x=rep(NA,2),pp.shft=rep(NA,2),mad=rep(NA,2)) ,
     D2O = data.frame(sp.idx=rep(NA,2),sp.x=rep(NA,2),pp.idx=rep(NA,2),pp.x=rep(NA,2),pp.shft=rep(NA,2),mad=rep(NA,2)) )
   for(c.mol in mols){rownames(out.l[[c.mol]]) <- poss}

  # Extrema of currrent profile [(H2O,D2O)-list of (max,min)-lists of (loc,idx)-lists]
   aux.extr.l <- lapply(mols , function(arg){extrema.locs(cbind(prof.in[,c('idx','x')],y=prof.in[,arg]),smoothen=TRUE)})
     names(aux.extr.l) <- mols
 # End of prelims
 #################################################################################################

  # Start/end SPs and PPs for (H2O/D2O) profile (SP = leftmost peak; PP = principal peak = peak nearest to pure profile right peak)
   for(c.mol in mols)
   { # iterate over H2O and D2O profiles

    # Load profile for current molecular type [(N.t))x(idx,x,H2O/D2O) df]
     c.prof <- prof.in[,c('idx','x',c.mol)]

    # Index of SP for current profile
     out.l[[c.mol]][ , 'sp.x'] <- min(aux.extr.l[[c.mol]]$max$x) # position
     out.l[[c.mol]][ , 'sp.idx'] <- aux.extr.l[[c.mol]]$max$idx[which.min(aux.extr.l[[c.mol]]$max$x)] # index

    for(c.pos in poss)
    { # iterate over start/end pure profiles

     # Index of pure PP (of current position and molecular type)
      aux.pure.idx <- pure.loc.l[[c.pos]][[c.mol]]['Rmax','idx']

     # Row hosting principal peak (PP) for current profile
      aux.row <- which.min(abs(aux.extr.l[[c.mol]]$max$idx - aux.pure.idx))
        if(aux.row < which.min(aux.extr.l[[c.mol]]$max$idx)){stop('ABORT: leftmost peak is principal peak [dat.vs.pure.R]')} # sanity check: PP =/= SP

     # Record PP index and position
      out.l[[c.mol]][c.pos,c('pp.idx','pp.x')] <- aux.extr.l[[c.mol]]$max[aux.row,c('idx','x') ]

     # Record PP shift wrt pure PP (of same position and molecular type)
      out.l[[c.mol]][c.pos,'pp.shft'] <- out.l[[c.mol]][c.pos,'pp.idx'] - aux.pure.idx # record principal peak shift

     # Crop pure profile [(idx,x,H2O/D2O) df]
      aux.idxs <- data.frame(min=pure.loc.l[[c.pos]][[c.mol]]['M','idx'] , max=pure.loc.l[[c.pos]][[c.mol]]['R','idx']) # pure index bounds [num vector]
      c.pure <- pure.prof.l[[c.pos]][pure.prof.l[[c.pos]]$idx %in% aux.idxs$min:aux.idxs$max , c('idx','x',c.mol)]
      c.pure[,c.mol] <- c.pure[,c.mol] - min(c.pure[,c.mol]) # remove constant background
      c.pure[,c.mol] <- c.pure[,c.mol]/(dt.g*sum(c.pure[,c.mol])) # renormalize cropped profile to unit AUC

     # Crop and shift input profile [(idx,x,H2O/D2O) df]
      c.dat <- prof.in[prof.in$idx %in% (aux.idxs$min:aux.idxs$max + out.l[[c.mol]][c.pos,'pp.shft']) , c('idx','x',c.mol)]
      c.dat[,c.mol] <- c.dat[,c.mol] - min(c.dat[,c.mol]) # remove constant background
      c.dat[,c.mol] <- c.dat[,c.mol]/(dt.g*sum(c.dat[,c.mol])) # normalize cropped profile to unit AUC

     # Mean absolute deviation between (cropped, renormalized) pure profile and data profile
      out.l[[c.mol]][c.pos,'mad'] <- sum(abs(c.dat[,c.mol] - c.pure[,c.mol]))/(aux.idxs$max-aux.idxs$min+1)

    } # end of iteration over start/end pure profiles
   } # end of iteration over H2O and D2O profiles

  if(length(unique(out.l[[c.mol]]$pp.idx))!=1){stop('ABORT: distinct start/end principal peaks [dat.vs.pure.R]')} # start/end principal peaks must match

    return(out.l)

} # end of function