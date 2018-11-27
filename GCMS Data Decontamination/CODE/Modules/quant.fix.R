####################################################################################################
# 'quant.fix.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This routine repairs (i.e. removes the pollutant contribution from) all datasets corresponding to
# the organ current in the global environment (set by settings.R).
#
# The removal is based on the decomposition f.dat(t) = prp*f.pure(t) + (1-prp)*f.poll(t) for the
# D2O dataset (profile), with the time t taking values from the rough interval grid and the scalar
# prp E [0,1] quantifying the (proportion of) "true signal" in the dataset. Here, f.dat, f.pure and
# f.poll are distributions derived from the dataset, pure profile and pollutant profile (for the
# current organ); as a result, they have unit AUC values. Although a (somehow) optimal choice of prp
# could be easily found (e.g. by brute force), 'smart cutoffs' create difficulties in practice.
# Specifically, the misalignment of the dataset's peaks with those of the pure/pollutant profiles
# complicates the decomposition. As a workaround, I avoid splitting a profile into signal/pollutant
# and only split, instead, its AUC over left/right intervals. Concretely, pure/pollutant profiles
# have their AUCs distinctly differently distibuted between left/right intervals: for pure profiles,
# most of AUC falls in the right interval; for pollutant profiles, nearly all of it is contained in
# the left one. As a result, any profile with an intermediate AUC splitting can be interpreted as an
# interpolation between pure/pollutant profiles. In detail, AUC is a linear functional of profiles,
# hence the decomposition above yields (for the AUCs over any sub-interval of the rough interval):
# AUC.dat = prp*AUC.pure + (1-prp)*AUC.poll. Both left and right sub-interval can be used to estimate
# prp; the result is the same (trivial to prove since AUC's sum to one).
#
# Given a dataset, all D2O profiles are first restricted (cropped) to that dataset's D2O interval
# (stored in the corresponding element of prof.loc.l) and renormalized to unit AUC. The splitting of
# AUC.poll (i.e. of one) into left/right AUC's (AUC.poll$L/AUC.poll$R) is trivial, given a pollutant
# profile. (If no pollutant profile's been tabulated, I assume that AUC.poll$L=1 and AUC.poll$R=0.)
# The best matching (among start/end) pure profile AUC is split similarly. As for the dataset, its
# unit AUC is split into AUC.dat$L and AUC.dat$R but in duplicate, as there are two D2O columns.
#
# Given the splittings AUC.poll$L and AUC.poll$R = 1 - AUC.poll$L (for the pollutant), AUC.pure$L
# and AUC.pure$R = 1 - AUC.pure$L (for the best-matching pure profile) and the duplicate splittings
# AUC.dat[typ,L] and AUC.poll[typ,R] = 1 - AUC.poll[typ,L] (with typ being D2O.A or D2O.B), I find
# prp (in duplicate) by solving AUC.dat$L = prp*AUC.pure$L + (1-prp)*AUC.poll$L. The solution
# is, naturally, prp = (AUC.dat$L - AUC.poll$L)/(AUC.pure$L - AUC.poll$L). These (duplicate) values
# are stored in the row of the df c.prp.df bearing the name of the dataset; its columns are named
# after the entries in cols$D2O.
#
# The D2O counts per interval are finally corrected and the corresponding (duplicate and average)
# ratios computed. Since we consider the signal to be prp*f.pure(t), the AUC of the dataset over an
# interval I will be AUC.dat.fix(I) = prp*AUC.pure(I). Naturally, the actual AUC of the profile over
# its entire 'smart' interval Io is not one but, instead, AUC.dat(Io). It follows that the AUC (over
# the interval I) of the dataset due to signal is AUC.dat.tot.fix(I) = prp*AUC.pure(I)*AUC.dat(Io),
# where prp E [0,1] the proportion of signal in the dataset, AUC.pure(I) E [0,1] the part of the
# (unit) AUC of the pure profile that corresponds to the interval I, and AUC.dat(Io) the AUC of the
# dataset over its entire sliding interval. With this formula, we repair the D2O counts per dataset
# (element of fils.v), duplicate (A,B) and interval (L,R,T). The ratios are then computed, for each
# dataset, per interval and duplicate, as is their average for each interval.
#
#   INPUTS
# None.
#
#   OUTPUTS
# * cnt.fix.l = filename-named list of (L,T,R)x(dat.cols,A,B,AB) df's; each df reports counts and
#                 D2O-over-H2O ratios for the dataset of the same name, summed up over the left (L),
#                 (L), right (R) and entire (T) interval; the four columns in dat.cols correspond to H2O and D2O
#                 counts in duplicate (see settings.R for order); 'A' and 'B' list the corresponding
#                 duplicate ratios and 'AB' the mean of these two ratios
#                                                [fils.v-list of (L,T,R)x(dat.cols,A,B,AB) num df's]
#
#   DEPENDENCIES
# <- settings.R
# <- presets.R
# <- all modules before it in main.R
# -> all modules after it in main.R
####################################################################################################
# Current (start/end) pure profiles [(start,end)-list of (N.t)x(idx,x,y) df's]
 c.pure <- pure.prof.l # current pure profiles [(start,end)-list of (N.t)x(idx,x,H2O,D2O) df's]
 c.pure <- lapply(c.pure , function(arg){arg[,-which(colnames(arg)=='H2O')]}) # remove H2O column [(start,end)-list of (N.t)x(idx,x,D2O) df's]
 for(c.pos in poss){colnames(c.pure[[c.pos]])[colnames(c.pure[[c.pos]])=='D2O'] <- 'y'} # rename D2O-col into y-col

# Remove pollutant individually from each (non-pollutant) dataset and record D2O (duplicate) counts and ratios into cnt.fix.l
 for(c.fil in dat.fils.df$fil)
 { # iterate over non-pollutant datasets

  # Load profiles for current dataset [(N.t)x(idx,x,dat.cols) df]
   c.cnt.df <- dat.tab.l[[c.fil]]

  # Cutoffs for current dataset and pure profile and BMMP positions
   c.LMR <- prof.loc.l[[c.fil]]$idx # dataset cutoffs [(L,M,R)x(dat.cols) df]
   c.pure.LMR <- pure.crop.l[[c.fil]]$idx # BMPP push-fwd cutoffs [(L,M,R)x(dat.cols) df]
   c.poss <- pure.crop.l[[c.fil]]$orig # BMPP positions (start or end) [(1)x(dat.cols) df]

  for(c.dupl in dupls)
  { # iterate over dataset duplicates

   # D2O data column for current duplicate (D2O.A or D2O.B) (string)
    c.dat.col <- intersect(cols[[c.dupl]],cols$D2O)

   # Crop pollutant profile to push-fwd D2O interval of current duplicate, remove constant offset and renormalize to unit AUC [(idx,x,y)-df]
    aux.poll <- subset(poll.prof.df , c.pure.LMR['L',c.dat.col] <= idx & idx <= c.pure.LMR['R',c.dat.col]) # crop current pollutant to push-fwd interval
    if(any(aux.poll$y>0)){aux.poll$y <- aux.poll$y/(dt.g*sum(aux.poll$y))} # renormalize cropped, offset-removed pollutant profile to unit AUC
   # Distribute AUC of pollutant profile to left/right intervals [1-D (L,R) df]
    c.AUC.poll <- data.frame( # left/right pollutant AUC's
      L = dt.g*sum(subset(aux.poll , idx < c.pure.LMR['M',c.dat.col])$y) ,
      R = dt.g*sum(subset(aux.poll , c.pure.LMR['M',c.dat.col] <= idx)$y) )
    if(all(c.AUC.poll==0)){c.AUC.poll <- data.frame(L=1 , R=0)} # left-localized pollutant, if there are no "pure pollutant" profiles

   # Crop BMPP to push-fwd D2O interval of current duplicate, remove constant offset and renormalize to unit AUC [(idx,x,y)-df]
    aux.pure <- subset(c.pure[[c.poss[[c.dat.col]]]] , c.pure.LMR['L',c.dat.col] <= idx & idx <= c.pure.LMR['R',c.dat.col])
    aux.pure$y <- aux.pure$y - min(aux.pure$y) # remove offset (computed over the cropped push-fwd interval of the curent D2O duplicate)
    aux.pure$y <- aux.pure$y/(dt.g*sum(aux.pure$y)) # renormalize cropped, offset-removed BMPP to unit AUC
   # Distribute AUC of BMPP to left/right intervals [1-D (L,R) df]
    c.AUC.pure <- data.frame( # left/right BMPP AUC's
      L = dt.g*sum(subset(aux.pure , idx < c.pure.LMR['M',c.dat.col])$y) ,
      R = dt.g*sum(subset(aux.pure , c.pure.LMR['M',c.dat.col] <= idx)$y) )

   # Crop current D2O duplicate to its interval, remove constant offset and renormalize to unit AUC [(idx,x,y)-df]
    aux.dat <- data.frame(idx=c.cnt.df$idx , x=c.cnt.df$x , y=c.cnt.df[,c.dat.col]) # remove irrelevant data cols [(idx,x,y)-df]
    aux.dat <- subset(aux.dat , c.LMR['L',c.dat.col] <= idx & idx <= c.LMR['R',c.dat.col])
    aux.dat$y <- aux.dat$y - min(aux.dat$y) # remove offset (computed over the interval of the curent D2O duplicate)
    aux.dat$y <- aux.dat$y/(dt.g*sum(aux.dat$y)) # renormalize cropped, offset-removed data profile to unit AUC
   # Distribute AUC of current D2O profile (in duplicate) to sliding left/right intervals [(cols$D2O)x(L,R) df]
    c.AUC.dat <- data.frame( # left/right data profile AUC's
      L = dt.g*sum(subset(aux.dat , idx < c.LMR['M',c.dat.col])$y) ,
      R = dt.g*sum(subset(aux.dat , c.LMR['M',c.dat.col] <= idx)$y) )

   # Record prp (per dataset; in duplicate) [(fils.v)x(cols$D2O) df]
    c.prp.df[c.fil,c.dat.col] <- min(max((c.AUC.dat$R - c.AUC.poll$R)/(c.AUC.pure$R - c.AUC.poll$R),0) , 1) # prp (restricted to [0,1])

  } # end of iteration over dataset duplicates

  # Repair D2O counts of current dataset: split prp-scaled profile AUC as in pure profile [(L,R,T)x(H2O.A,D2O.A,H2O.B,D2O.B,A,B,AB) df]
   c.cnt.fix <- cnt.l[[c.fil]] # initialize with unrepaired counts (H2O counts will not change below)
   for(c.intv in ratio.rows)
   { # loop over left/right/entire interval
    for(c.dat.col in cols$D2O)
    { # repair counts
     c.cnt.fix[c.intv,c.dat.col] <- round(
       ifelse(c.intv=='T',1,c.AUC.pure[[c.intv]])* # pure profile AUC for current interval (one, for entire interval)
       c.prp.df[c.fil,c.dat.col]* # prp-value for current D2O duplicate
       cnt.l[[c.fil]][['T',c.dat.col]] ) # total profile AUC (i.e. over entire interval) for current duplicate
    }
    c.cnt.fix[c.intv,'A'] <- c.cnt.fix[c.intv,cols$A[cols$A %in% cols$D2O]]/c.cnt.fix[c.intv,cols$A[cols$A %in% cols$H2O]] # recompute duplicate A ratio
    c.cnt.fix[c.intv,'B'] <- c.cnt.fix[c.intv,cols$B[cols$B %in% cols$D2O]]/c.cnt.fix[c.intv,cols$B[cols$B %in% cols$H2O]] # recompute duplicate B ratio
    c.cnt.fix[c.intv,'AB'] <- rowMeans(c.cnt.fix[c.intv,cols$D2O]/c.cnt.fix[c.intv,cols$H2O]) # recompute duplicate-averaged ratio
   }

  # Update corresponding cnt.fix.l entry with (and, if write switch is ON, also write to disk) df holding repaired counts
   cnt.fix.l[[c.fil]] <- c.cnt.fix # assign repaired counts df to element of cnt.fix.l bearing current dataset's name
   if(write.swtc){write.csv(c.cnt.fix , paste('../Results/',organ,'/CSV_Results/',c.fil,sep=''))} # only (re)write outcome if write switch is ON
     
 } # end of iterate over non-pollutant datasets

# Bind signal proportions of datasets into SR.df, then (if write switch is ON) write SR.df to disk
 SR.df <- setNames(cbind(SR.df,c.prp.df) , c(colnames(SR.df),dupls))
 if(write.swtc){write.csv(SR.df , paste('../Results/',organ,'/CSV_Results/SR.csv',sep=''))} # only (re)write outcome if write switch is ON