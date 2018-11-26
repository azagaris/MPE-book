####################################################################################################
# 'quant.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This routine returns AUCs for all datasets corresponding to the organ current in the workspace.
# AUCs are computed over three distinct 'smart cutoff intervals': the left interval, ranging from
# a (dataset-specific) left cutoff to the valley, the right interval (valley to right cutoff) and
# the entire ('total') interval (left-to-right cutoff).
#
# Computing AUCs is a mattter of truncating datasets into their smart cutoff intervals, removing (or
# not) a constant offset (i.e. the dataset minimum over the smart interval) and summing up counts.
# (I do not multiply by dt.g, since counts are not rates; multiplying or not does not affect the
# computed deuteration ratios.) Smart cutoffs are retrieved (per dataset) from the list prof.loc.l,
# cf. dat.locs.R. Datasets are retrieved from the list dat.tab.l, cf. presets.R. The cardinalities
# of the truncated H2O/D2O datasets differ, since the H2O and D2O cutoffs vary per dataset and even
# duplicate, and hence cropped H2O/D2O datasets cannot be stored side by side as columns of some df.
# Instead, I store them in the list c.cnt.l, having elements (H2O,D2O) which are themselves df's
# with columns (idx,x,H2O.A,H2O.B) and (idx,x,D2O.A,D2O.B), respectively. This list is rewritten as
# we iterate over datasets within the current organ.
#
# T(otal)-counts (per molecular type and duplicate) are obtained by summing up data columns in each
# element of c.cnt.l. These are stored in the T-row of the dataset-specific df c.cnt, in the columns
# named after dat.cols (cf. settings.R). Also stored in the same row are deuteration ratios, namely
# in columns 'A' and 'B' (for count(D2O.A)/count(H2O.A) and count(D2O.B)/count(H2O.B)), as well as
# their average (column 'AB'). The process is repeated with left and right half-intervals and the
# results stored in the 'L' and 'R' rows of c.cnt. When ready, the df is stored as the homonymous
# element of the list cnt.l.
#
# Next to AUCs, we also compute (per dataset) a 'skewness ratio' (SR) estimating the presence of
# pollutant (in each dataset). For datasets in the control ensemble, where pollutant is absent, the
# proportion of AUC contained in the left interval is rather stable, i.e. its variability is
# markedly smaller than for non-control datasets. Since contaminant is localized nearly exclusively
# in the left interval, increased pollutant presence "skews" that ratio towards greater values. That
# SR is computed per dataset in a series of steps: first, we compute the left-over-total-interval
# AUC ratio; second, we subtract from it a reference SR-value (equaling the mean SR over all control
# datasets); and finally, we express the result in sd-units, where 'sd' is the standard deviation of
# SR's across control datasets. In this setting, a negative (positive) skewness ratio signifies a
# left-interval-AUC-proportion below (above) the reference ratio. In practice, the occasional D2O
# dataset shows a moderate negative value of a few sd, but most datasets exhibit positive values in
# the tens of sd; H2O datasets only show moderate positive and negative skewness ratio values.
#
#   INPUTS
# None.
#
#   OUTPUTS
# * cnt.l = filename-named list of (L,R,T)x(dat.cols,A,B,AB) df's; each df reports the counts and
#           D2O-over-H2O ratios for the dataset of the same name, summed up over the left (L), right
#           (R) or entire (T) interval; the four dat.cols-columns correspond to H2O and D2O counts
#           in duplicate (see settings.R for order), whereas A and B list the ratios per duplicate
#           and AB the means of these two
#                                                [fils.v-list of (L,R,T)x(dat.cols,A,B,AB) num df's]
#
# * SR.df = finename-named df whose rows are named after the datasets corresponding to the current
#           organ and whose columns are dat.cols-named; each entry is a skewness ratio (see above
#           for details) for the homonymous dataset column [(fils.v)x(dat.cols) df's]
#
#   DEPENDENCIES
# <- settings.R
# <- presets.R
# <- all modules before it in main.R
# -> all modules after it in main.R
####################################################################################################
# Dataset quantification
 for(c.fil in dat.fils.df$fil)
 { # iterate over non-pure-pollutant datasets for current organ

  # Load profiles for current dataset [(N.t)x(idx,x,dat.cols) df]
   c.cnt.df <- dat.tab.l[[c.fil]]

  # Load left/mid/right cutoff indices per duplicate of current dataset [(L,M,R)x(dat.cols) df]
   c.LMR <- prof.loc.l[[c.fil]]$idx

  # Make (H2O,D2O)-list to hold trimmed current dataset [dat.cols-list of (idx,cnt)-col df's]
   c.cnt.l <- setNames( vector('list',length(dat.cols)) , dat.cols)

  # Trim each molecular type to fit left/right cutoffs (specific to dataset), remove offset (or not) and store in c.cnt.l-element
   for(c.dat.col in dat.cols)
   { # iterate over data columns (i.e. over molecular types and duplicates)
    c.idx <- data.frame(t(c.LMR[,c.dat.col,drop=F])) # L/M/R cutoffs (indices) for current data col [(L,M,R)-col df]
    c.cnt.l[[c.dat.col]] <- c.cnt.df[c.idx$L <= c.cnt.df$idx & c.cnt.df$idx <= c.idx$R , c('idx',c.dat.col)]
    c.cnt.l[[c.dat.col]][,c.dat.col] <- c.cnt.l[[c.dat.col]][,c.dat.col] - min(c.cnt.l[[c.dat.col]][,c.dat.col]) # remove offset
   } # end of iteration over data columns

  # Record totql counts over L/R/T intervals in c.cnt [(L,R,T)x(dat.cols,dupls,AB) df]
   c.cnt <- setNames(data.frame(matrix(NA,length(ratio.rows),length(ratio.cols)) , row.names=ratio.rows) , ratio.cols) # NA-preallocate counts df [(L,R,T)x(dat.cols) df]
   for(c.dat.col in dat.cols)
   { # compute counts per data column (i.e. per molecular types and duplicates)
    c.idx <- data.frame(t(c.LMR[,c.dat.col,drop=F])) # L/M/R cutoffs (indices) for current data col [(L,M,R)-col df]
    aux.cnt.df <- c.cnt.l[[c.dat.col]] # (cropped) count df for current data col [(idx,dat.col)-col df]
    c.cnt['L',c.dat.col] <- sum(aux.cnt.df[aux.cnt.df$idx < c.idx$M , c.dat.col]) # L-interval AUC
    c.cnt['R',c.dat.col] <- sum(aux.cnt.df[c.idx$M <= aux.cnt.df$idx , c.dat.col]) # R-interval AUC
    c.cnt['T',c.dat.col] <- sum(aux.cnt.df[,c.dat.col]) # T-interval AUC (sum of L- and R-AUCs)
   }

  # Record deuteration ratios over L/R/T intervals in c.cnt [(L,R,T)x(dat.cols,dupls,AB) df]
   for(c.row in ratio.rows)
   { # D2O/H2O count ratios (in duplicate) per interval
    for(c.dupl in dupls) # deuteration ratio per duplicate
      {c.cnt[c.row,c.dupl] <- c.cnt[c.row , cols[[c.dupl]][cols[[c.dupl]] %in% cols$D2O]]/c.cnt[c.row , cols[[c.dupl]][cols[[c.dupl]] %in% cols$H2O]]}
    c.cnt[c.row,ratio.cols[!ratio.cols %in% c(dat.cols,dupls)]] <- rowMeans(c.cnt[c.row,dupls]) # mean ratio
   }

  # Update corresponding cnt.l entry
   cnt.l[[c.fil]] <- c.cnt # assign counts df to element of cnt.l bearing current dataset's name

  # Record H2O and D2O (in duplicate) skewness ratios
   SR.df[c.fil,] <- apply(c.cnt[,dat.cols] , 2 , function(arg){arg['L']/arg['T']})

 } # end of iteration over non-pure-pollutant datasets for current organ

# Skewness ratios statistics for control ensembles [(mean,sd)x(H2O,D2O) df]
 aux.SR.ctrl.stats <- data.frame( # statistics df
   H2O = c(mean=mean(as.matrix(SR.df[ctrl.fils.df$fil , cols$H2O])) , sd=sd(as.matrix(SR.df[ctrl.fils.df$fil , cols$H2O]))) ,
   D2O = c(mean=mean(as.matrix(SR.df[ctrl.fils.df$fil , cols$D2O])) , sd=sd(as.matrix(SR.df[ctrl.fils.df$fil , cols$D2O]))) )

# Recast dataset skewness ratios as aberrations (from the control-ensemble mean) measured in sd-units [(c.fils.v)x(dat.cols) df]
 SR.df <- sweep( sweep(SR.df,2,unlist(aux.SR.ctrl.stats['mean',]),'-') , 2 , unlist(aux.SR.ctrl.stats['sd',]) , '/' )