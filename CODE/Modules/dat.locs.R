####################################################################################################
# 'dat.locs.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This routine determines 'smart (cutoff) intervals' for all datasets (save for pure-pollutant ones)
# within the current organ. All datasets were initially cropped to an organ-specific rough interval
# (set in settings.R). For an accurate analysis, this must be cropped further and per profile, i.e.
# independently for each duplicate of each molecular type (undeuterated/deuterated). This is done
# to avoid counting irrelevant particles hitting the detector at nearby retention times.
#
# These smart intervals are based on the cutoff intervals for the (start/end) pure profiles of the
# current organ, which contain 100*(1-peak.thr)% of the pure profile AUC; see pure.prof.locs.R.
# Since dataset (i.e. neither pure nor pure-pollutant) profiles are generally polluted, it is not
# necessarily wise to work with the AUC as done for the pure profiles. Furthermore, dataset profiles
# contain impurities, manifesting themselves as (low but substantial) peaks near signal peaks; our
# smart intervals should exclude such impurities to the extent possible.
#
# To that end, and given a dataset (the module iterates over all datasets), we align the start/end
# pure profiles with the given dataset profile, i.e. shift them horizontally until their principal
# peaks (PPs) align with the given profile's PP. (Each dataset is composed of an undeuterated (H2O)
# and a deuterated (D2O) profile - in duplicate, too; thus, alignment etc is performed in duplicate,
# for both H2O/D2O profiles (using the H2O/D2O pure profiles) and for both start/end pure profiles.
# Also, pure profile PP means, here, the right peak holding most of the AUC; the definition of a PP
# for a dataset profile is more involved, see the routine dat.vs.pure.R performing the alignment.)
# The shifted left (L)/right (R) cutoffs of the (start/end) pure profiles yield 'cutoff candidates'
# for the L/R cutoffs of the dataset profile. To choose between start/end-derived candidates, we
# compute the mean average deviations (MADs) between profile and aligned pure profiles, then select
# the pure profile yielding the smallest MAD; again, see dat.vs.pure.R for details.
#
# For each duplicate, these L/R candidates occupy two entries in the 4x(mol,orig,idx,x,dist) df's
# c.cand.L and c.cand.R, respectively. In these entries, 'mol' is H2O (in one)/D2O (in the other)
# (marking whether the cutoff applies to H2O/D2O profile of the duplicate), 'orig' is either 'start'
# or 'end' (marking the best match, which also yields the cutoffs), 'idx' and 'x' encode the index-
# and x-values of the cutoff (as a gridpoint), and 'dist' is a natural number measuring the distance
# (in dt.g units) of the cutoff candidate from either PP (in c.cand.R) or secondary peak (SP) (in
# c.cand.L).
#
# The other pair of cutoff candidates derives from local minima ('vallleys'). For L-cutoffs, we find
# all minima below SP and retain the rightmmost one; for R-cutoffs, we retain the leftmost among all
# minima above PP. These are recorded in the remaining two rows of c.cand.L and c.cand.R, with
# 'orig' set to 'valley' and all other columns having the same meaning as before.
#
# The actual L-/R-cutoffs are derived from these candidates; we cover R-cutoffs here, remarking that
# L-cutoffs are derived analogously. We call a profile pure-limited if the minima ('valleys') above
# its PP (if any) also fall above its pure-derived R-cutoff (obtained by shifting the closest
# matching pure profile so that its PP aligns with the dataset PP). We call it valley-limited if the
# leftmost valley above its PP falls below (i.e. closer to the PP) than its pure-derived R-cutoff.
# Now, if both H2O/D2O profiles (of the curent duplicate) are pure-limited, then the pure-derived
# candidates become R-cutoffs. If one of the profiles is pure-limited but the other valley-limited,
# the R-cutoff of the latter is set to the corresponding valley. As for the former, its pure-derived
# candidate is shifted closer to the PP and set as cutoff, to mimic the tail excised from the latter
# profile. Specifically, if the distance (for the latter) from PP to valley candidate is 0<A<1 of
# that from PP to pure-derived candidate, then the pure-derived candidate for the former is shifted
# so that the ratio of distances from PP to new and old pure-derived cutoff is also A. Finally, if
# both profiles are valley-limited, then we follow the same procedure but upon selecting the profile
# that is most stringently limited (i.e. has the smallest A, in the above notation).
#
# The routine also computes the midpoint of each profile, i.e. the point separating the two modes
# (left/right). This is found as the rightmost minimum between SP and PP, and we select specifically
# the rightmost because minima to its left might be due to bimodality created by interfrence of the
# pollutant with the SP.
#
#   INPUTS
# None.
#
#   OUTPUTS
# * pure.crop.l = ...
#
# * prof.loc.l = dat.fils.v-named list, each element of which corresponds to the eponymous dataset
#                and is itself a (idx,x)-list; the idx (x) element of that list records the index
#                (x-value) of the cutoffs (computed as detailed above) for all four profiles in the
#                dataset (H2O/D2O in duplicate); the cutoffs themselves are labeled L, M and R for
#                left, midpoint and right
#                                 [dat.fils.v-list of (idx,x,orig)-lists of (L,M,R)x(dat.cols) df's]
#
#   DEPENDENCIES
# <- settings.R
# <- presets.R
# <- pure.profs.tab.R
# <- extrema.locs.R
# -> all modules after it in main.R
####################################################################################################
for(c.fil in dat.fils.df$fil)
{ # iterate over datasets for current organ

 # Load profiles for current dataset (counts over regular grid on rough interval) [(N.t)x(idx,x,dat.cols) df]
  c.cnt.df <- dat.tab.l[[c.fil]]

 # Iterate over current dataset duplicates
  for(c.dupl in dupls)
  { # iterate over dataset duplicates

   # Preliminaries
    # Preallocate df's to hold right and left cutoff candidates for H2O/D2O profiles [(mol,orig,idx,x,dist) df's]
     c.cand.R <- data.frame(mol=c(rep(mols[1],2),rep(mols[2],2)) , orig=rep(c('valley','pure'),2) , idx=rep(Inf,4) , x=rep(Inf,4) , dist=rep(Inf,4))
         levels(c.cand.R$orig) <- levels(c.cand.R$orig) <- c(levels(c.cand.R$orig),poss) # additional levels replace 'pure' below
     c.cand.L <- c.cand.R # left cutoff candidate df
       c.cand.L[,c('idx','x')] <- -Inf # Inf-entries -> -Inf
    # Store (H2O/D2O) count-profiles for current duplicate [(N.t)x(idx,x,H2O,D2O) df]
     aux.df <- c.cnt.df[ , c('idx','x',cols[[c.dupl]])] # excise non-current duplicate
       for(c.col in cols[[c.dupl]]){names(aux.df)[names(aux.df)==c.col] <- strsplit(c.col,'[.]')[[1]][1]} # rename data cols
    # Store extrema of (H2O/D2O) count-profiles for current duplicate [(H2O,D2O)-list of (max,min)-lists of (loc,idx)-lists]
     c.extr.l <- lapply(mols , function(arg){extrema.locs(cbind(aux.df[,c('idx','x')],y=aux.df[,arg]),smoothen=TRUE)})
       names(c.extr.l) <- mols
    # Comparison between current duplicate and start/end pure profiles [mols-list of (start,end)x(idx,x,shft,mad) df's]
     c.peak.l <- dat.vs.pure(aux.df)
    # Record position (i.e. start/end) of best-matching pure profiles (BMPPs) for H2O/D2O datasets of current duplicate
     for(c.mol in mols)
      {pure.crop.l[[c.fil]]$orig[[intersect(cols[[c.mol]],cols[[c.dupl]])]] <- rownames(c.peak.l[[c.mol]])[which.min(c.peak.l[[c.mol]]$mad)]}

   # Midpoint: find minima between leftmost and primary peaks, then set midpoint to rightmost among them
    for(c.mol in mols)
    { # iterate over H2O and D2O profiles for current duplicate
     c.dat.col <- intersect(cols[[c.mol]],cols[[c.dupl]]) # dataset corresponding to current molecular type and duplicate [dat.cols-element]
     c.v <- min(c.extr.l[[c.mol]]$max$x)<c.extr.l[[c.mol]]$min$x & c.extr.l[[c.mol]]$min$x<unique(c.peak.l[[c.mol]]$pp.x) # 'is min interpeak?' [logi vec]
     prof.loc.l[[c.fil]]$idx['M',c.dat.col] <- tail(c.extr.l[[c.mol]]$min$idx[c.v] , 1) # store index
     prof.loc.l[[c.fil]]$x['M',c.dat.col] <- tail(c.extr.l[[c.mol]]$min$x[c.v] , 1) # store location
     c.pos <- pure.crop.l[[c.fil]]$orig[[c.dat.col]] # BMPP position [start or end]
     pure.crop.l[[c.fil]]$idx['M',c.dat.col] <- pure.loc.l[[c.pos]][[c.mol]]['M','idx'] # BMPP midpoint (index)
     pure.crop.l[[c.fil]]$x['M',c.dat.col] <- pure.loc.l[[c.pos]][[c.mol]]['M','x'] # BMPP midpoint (x-value)
    } # end of iteration over H2O and D2O profiles for current duplicate

   # Right cutoff candidates
    for(c.mol in mols)
    { # iterate over molecular types (H2O/D2O) of current duplicate
     # Indices (in c.extr.l[[c.mol]]$min$x/idx) of profile min above profile PP [nat vector]
      c.v <- which(c.extr.l[[c.mol]]$min$x > unique(c.peak.l[[c.mol]]$pp.x))
     # If profile has minima above its PP, then record leftmost such min (valley) as right cutoff candidate
      if(length(c.v)>=1)
      { # check whether any minima are above PP
       c.x <- min(c.extr.l[[c.mol]]$min$x[c.v]) # location of leftmost min above PP [scalar]
       c.idx <- c.extr.l[[c.mol]]$min$idx[ c.v[which.min(c.extr.l[[c.mol]]$min$x[c.v])] ] # index of leftmost min above PP [nat]
       c.dist <- c.idx - unique(c.peak.l[[c.mol]]$pp.idx) # distance (from profile PP) of leftmost min above PP (in dt.g units) [nat]
       c.cand.R[c.cand.R$mol==c.mol & c.cand.R$orig=='valley' , c('idx','x','dist')] <- data.frame(idx=c.idx,x=c.x,dist=c.dist)
      } # end of check whether any minima are above PP
     # Also record (shifted) right cutoff of best-matching pure profile as right cutoff candidate
      c.pos <- rownames(c.peak.l[[c.mol]])[which.min(c.peak.l[[c.mol]]$mad)] # best-matching profile ['start' or 'end']
      c.idx <- pure.loc.l[[c.pos]][[c.mol]]['R','idx'] + c.peak.l[[c.mol]][c.pos,'pp.shft'] # index of shifted, pure-derived right cutoff
      c.x <- c.cnt.df[c.cnt.df$idx==c.idx,'x'] # location of shifted, pure-derived right cutoff
      c.dist <- c.idx - unique(c.peak.l[[c.mol]]$pp.idx) # distance (from profile PP) of shifted, pure-derived right cutoff (in dt.g units)
      c.cand.R[c.cand.R$mol==c.mol & c.cand.R$orig=='pure' , c('orig','idx','x','dist')] <- data.frame(orig=c.pos,idx=c.idx,x=c.x,dist=c.dist)
    } # end of iterate over molecular types (H2O/D2O) of current duplicate

   # Left cutoff candidates
    for(c.mol in mols)
    { # iterate over molecular types (H2O/D2O) of current duplicate
     # Indices (in c.extr.l[[c.mol]]$min$x/idx) of profile min below profile SP [nat vector]
      c.v <- which(c.extr.l[[c.mol]]$min$x < min(c.extr.l[[c.mol]]$max$x))
     # If profile has minima below its SP, then record rightmost such min (valley) as left cutoff candidate
      if(length(c.v)>=1)
      { # check whether any minima are below SP
       c.x <- max(c.extr.l[[c.mol]]$min$x[c.v]) # location of rightmost min below SP [scalar]
       c.idx <- c.extr.l[[c.mol]]$min$idx[ c.v[which.max(c.extr.l[[c.mol]]$min$x[c.v])] ] # index of rightmost min below SP [nat]
       c.dist <- unique(c.peak.l[[c.mol]]$sp.idx) - c.idx # distance (from profile SP) of rightmost min below SP (in dt.g units) [int]
       c.cand.L[c.cand.L$mol==c.mol & c.cand.L$orig=='valley' , c('idx','x','dist')] <- data.frame(idx=c.idx,x=c.x,dist=c.dist)
      } # end of check whether any minima are below SP
     # Also record (shifted) left cutoff of best-matching pure profile as left cutoff candidate
      c.pos <- rownames(c.peak.l[[c.mol]])[which.min(c.peak.l[[c.mol]]$mad)] # best-matching profile ['start' or 'end']
      c.idx <- pure.loc.l[[c.pos]][[c.mol]]['L','idx'] + c.peak.l[[c.mol]][c.pos,'pp.shft'] # index of shifted, pure-derived left cutoff
      c.x <- c.cnt.df[c.cnt.df$idx==c.idx,'x'] # location of shifted, pure-derived left cutoff
      c.dist <- unique(c.peak.l[[c.mol]]$sp.idx) - c.idx # distance (from profile SP) of shifted, pure-derived left cutoff (in dt.g units)
      c.cand.L[c.cand.L$mol==c.mol & c.cand.L$orig=='pure' , c('orig','idx','x','dist')] <- data.frame(orig=c.pos,idx=c.idx,x=c.x,dist=c.dist)
    } # end of iterate over molecular types (H2O/D2O) of current duplicate

   # Record left/mid/right cutoff candidates
    for(c.mol in mols)
    { # iterate over molecular types (H2O/D2O) of current duplicate
      c.dat.col <- intersect(cols[[c.mol]],cols[[c.dupl]]) # dat.cols-element corresponding to current molecular type and duplicate [string]
     # Pure-profile right candidates
      prof.pure.cand.l[[c.fil]]$idx['R',c.dat.col] <- c.cand.R[c.cand.R$mol==c.mol & c.cand.R$orig!='valley','idx']
      prof.pure.cand.l[[c.fil]]$x['R',c.dat.col] <- c.cand.R[c.cand.R$mol==c.mol & c.cand.R$orig!='valley','x']
     # Super-PP valley right candidates
      prof.valley.cand.l[[c.fil]]$idx['R',c.dat.col] <- c.cand.R[c.cand.R$mol==c.mol & c.cand.R$orig=='valley','idx']
      prof.valley.cand.l[[c.fil]]$x['R',c.dat.col] <- c.cand.R[c.cand.R$mol==c.mol & c.cand.R$orig=='valley','x']
     # Pure-profile left candidates
      prof.pure.cand.l[[c.fil]]$idx['L',c.dat.col] <- c.cand.L[c.cand.L$mol==c.mol & c.cand.L$orig!='valley','idx']
      prof.pure.cand.l[[c.fil]]$x['L',c.dat.col] <- c.cand.L[c.cand.L$mol==c.mol & c.cand.L$orig!='valley','x']
     # Super-PP valley left candidates
      prof.valley.cand.l[[c.fil]]$idx['L',c.dat.col] <- c.cand.L[c.cand.L$mol==c.mol & c.cand.L$orig=='valley','idx']
      prof.valley.cand.l[[c.fil]]$x['L',c.dat.col] <- c.cand.L[c.cand.L$mol==c.mol & c.cand.L$orig=='valley','x']
    } # end of iterate over molecular types (H2O/D2O) of current duplicate

   # Determine right/left cutoffs for H2O/D2O profiles
    for(c.sid in c('L','R'))
    { # treat left/right cutoffs separately

     # Pass the current-side candidate df (c.cand.L or c.cand.R) to the new variable c.cand [(mol,orig,idx,x,dist) df]
      c.cand <- eval(parse(text=paste('c.cand.',c.sid,sep='')))
     # Initialize the current-side cutoff df [(mol,orig,idx,x,dist) df]
     c.cutoff <- data.frame(mol=mols , orig=c('pure','pure') , idx=c(NA,NA) , x=c(NA,NA) , dist=c(NA,NA)) # orig-entry initializes levels
       levels(c.cutoff$orig) <- levels(c.cand.R$orig) # fix levels
     # Select, for each (H2O/D2O) profile, its cutoff candidate closest to the relevant peak [(mol,orig,idx,x,shft) df]
      for(c.mol in mols){c.cutoff[c.cutoff$mol==c.mol,] <- c.cand[c.cand$mol==c.mol,][which.min(c.cand[c.cand$mol==c.mol,'dist']),]}

     # Case: if both H2O/D2O selected candidate cutoffs are pure-derived (no valleys), set them as cutoffs
      if(all(c.cutoff$orig!='valley'))
      {
       for(c.mol in mols)
       {
        prof.loc.l[[c.fil]]$idx[c.sid,intersect(cols[[c.mol]],cols[[c.dupl]])] <- c.cutoff[c.cutoff$mol==c.mol,'idx'] # store index
        prof.loc.l[[c.fil]]$x[c.sid,intersect(cols[[c.mol]],cols[[c.dupl]])] <- c.cutoff[c.cutoff$mol==c.mol,'x'] # store location
       }
      } # end of case: no profile exhibits valley

     # Case: if a single H2O/D2O selected candidate cutoff is a valley, set it as cutoff and shrink the other profile's cutoff accordingly
      if(any(c.cutoff$orig!='valley') & any(c.cutoff$orig=='valley'))
      {
       c.mol <- as.character(c.cutoff[c.cutoff$orig=='valley','mol']) # valley-ltd profile ['H2O' or 'D2O']
       prof.loc.l[[c.fil]]$idx[c.sid,intersect(cols[[c.mol]],cols[[c.dupl]])] <- c.cutoff[c.cutoff$mol==c.mol,'idx'] # set c.mol cutoff to valley (index)
       prof.loc.l[[c.fil]]$x[c.sid,intersect(cols[[c.mol]],cols[[c.dupl]])] <- c.cutoff[c.cutoff$mol==c.mol,'x'] # set c.mol cutoff to valley (location)
       aux.mol <- setdiff(mols,c.mol) # non-valley-ltd profile ['H2O' or 'D2O']
       c.rat <- subset(c.cand,mol==c.mol & orig=='valley')$dist/subset(c.cand,mol==c.mol & orig!='valley')$dist # shrink factor [scalar]
       aux.mol.dist.srnk <- round((1-c.rat)*c.cutoff[c.cutoff$mol==aux.mol,'dist']) # shrinkage of cutoff-to-peak distance for non-valley-ltd profile
       c.cutoff[c.cutoff$mol==aux.mol,'dist'] <- c.cutoff[c.cutoff$mol==aux.mol,'dist'] - aux.mol.dist.srnk # shrunk distance
       c.cutoff[c.cutoff$mol==aux.mol,'idx'] <- c.cutoff[c.cutoff$mol==aux.mol,'idx'] + ifelse(c.sid=='L',1,-1)*aux.mol.dist.srnk # new index
       c.cutoff[c.cutoff$mol==aux.mol,'x'] <- c.cnt.df[c.cnt.df$idx==c.cutoff[c.cutoff$mol==aux.mol,'idx'] , 'x'] # new location
       prof.loc.l[[c.fil]]$idx[c.sid,intersect(cols[[aux.mol]],cols[[c.dupl]])] <- c.cutoff[c.cutoff$mol==aux.mol,'idx'] # new cutoff (index)
       prof.loc.l[[c.fil]]$x[c.sid,intersect(cols[[aux.mol]],cols[[c.dupl]])] <- c.cutoff[c.cutoff$mol==aux.mol,'x'] # new cutoff (location)
      } # end of case: a single profile exhibits valley

     # Case: if both H2O/D2O selected candidate cutoffs are valleys, set shortest as cutoff and shrink the other profile's cutoff accordingly
      if(all(c.cutoff$orig=='valley'))
      {
       c.rats <- c(NA,NA)
         names(c.rats) <- mols
       for(c.mol in mols){c.rats[c.mol] <- subset(c.cand,mol==c.mol & orig=='valley')$dist/subset(c.cand,mol==c.mol & orig!='valley')$dist} # shrink factors
       c.mol <- names(which.min(c.rats)) # shortest valley-ltd profile ['H2O' or 'D2O']
       prof.loc.l[[c.fil]]$idx[c.sid,intersect(cols[[c.mol]],cols[[c.dupl]])] <- c.cutoff[c.cutoff$mol==c.mol,'idx'] # set c.mol cutoff to valley (index)
       prof.loc.l[[c.fil]]$x[c.sid,intersect(cols[[c.mol]],cols[[c.dupl]])] <- c.cutoff[c.cutoff$mol==c.mol,'x'] # set c.mol cutoff to valley (location)
       aux.mol <- setdiff(mols,c.mol) # other valley-ltd profile ['H2O' or 'D2O']
       c.rat <- c.rats[c.mol] # shrink factor [scalar]
       aux.mol.dist.srnk <- round((1-c.rat)*c.cand[c.cand$mol==aux.mol & c.cand$orig!='valley','dist']) # shrinkage of cutoff-to-peak distance for aux.mol
       c.cutoff[c.cutoff$mol==aux.mol,'dist'] <- subset(c.cand,mol==aux.mol & orig!='valley')$dist - aux.mol.dist.srnk # shrunk distance
       c.cutoff[c.cutoff$mol==aux.mol,'idx'] <- subset(c.cand,mol==aux.mol & orig!='valley')$idx + ifelse(c.sid=='L',1,-1)*aux.mol.dist.srnk # new index
       c.cutoff[c.cutoff$mol==aux.mol,'x'] <- c.cnt.df[c.cnt.df$idx==c.cutoff[c.cutoff$mol==aux.mol,'idx'] , 'x'] # new location
       prof.loc.l[[c.fil]]$idx[c.sid,intersect(cols[[aux.mol]],cols[[c.dupl]])] <- c.cutoff[c.cutoff$mol==aux.mol,'idx'] # new cutoff (index)
       prof.loc.l[[c.fil]]$x[c.sid,intersect(cols[[aux.mol]],cols[[c.dupl]])] <- c.cutoff[c.cutoff$mol==aux.mol,'x'] # new cutoff (location)
      } # end of case: both profiles exhibit valleys

   # Push forward and record right/left cutoffs for H2O/D2O BMPPs of current duplicate
    for(c.mol in mols)
    {
     c.dat.col <- intersect(cols[[c.mol]],cols[[c.dupl]]) # dataset corresponding to current molecular type and duplicate [dat.cols-element]
     c.pos <- pure.crop.l[[c.fil]]$orig[[c.dat.col]] # BMPP position [start or end]
     pure.crop.l[[c.fil]]$idx[c.sid,c.dat.col] <- c.cutoff[c.cutoff$mol==c.mol,'idx'] - c.peak.l[[c.mol]][c.pos,'pp.shft'] # push index fwd
    }

    } # end of treat left/right cutoffs separately

  } # end of iteration over dataset duplicates

} # end of iteration over datasets for current organ

    print(sprintf('Smart intervals for data profiles complete [%s]',organ))