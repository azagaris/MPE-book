####################################################################################################
# 'poll.tab.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This module tabulates the pure-pollutant profile for the organ set in settings.R. Contrary to pure
# profiles (derived from control profiles), there are no "pollutant controls" to estimate pollutant
# from. Instead, we choose datasets (see settings.R) in which the pollutant seems disparately large
# and obtain the pollutant by averaging/mollifying these. Importantly, the pure-pollutant profile is
# a D2O profiles there is little evidence of H2O profile contamination, hence we neither tabulate an
# H2O pollutant nor correct H2O datasets for contamination.
#
# The chosen pure-pollutant profiles are read from the fils.poll.l element (cf. settings.R) bearing
# the name of the organ currently in the workspace. That list element consists of the .csv filenames
# that will inform the construction of our pure-pollutant profile. The procedure itself is similar
# to that in pure.tabs.R: we first form the df c.poll.df, having an index column ('idx'), time grid
# column ('x') and twice (to accommodate duplicate measurements) as many extra columns as there are
# pure-pollutant datasets. These extra columns originally store the D2O counts of these datasets and
# are subsequently normalized so that their AUC equals one; this makes each such column a tabulated
# distribution. The average of these distributions is also a distribution (the average pollutant
# profile), which we subsequently mollify with a Gaussian kernel (see presets.R). This mollified
# distribution is the desired pollutant profile and saved in the df poll.prof.df.
#
#   INPUTS
# None.
#
#   OUTPUTS
# * poll.prof.df = df holding an index column, a timegrid column and a column tabulating the
#                  (average, mollified, deuterated) pollutant over the rough interval
#                                                                               [(N.t)x(idx,x,y) df]
#
#   DEPENDENCIES
# <- settings.R
# <- presets.R
# -> quant.fix.R
####################################################################################################
# List (in duplicate) normalized D2O profiles of datasets in fils.poll.v (skipped if fils.poll.v = NULL) [N.t-long df]
 c.poll.df <- data.frame(idx=idx.g , x=t.g) # initialize df holding pure-pollutant profile [(N.t)x(idx,x) df]
 for(c.fil in fils.poll.l[[organ]])
 { # iterate over "pure pollutant" datasets
  # Load D2O data (in duplicate) of the current dataset over rough interval [(N.t)x(cols$D2O) df]
   c.dat.df <- dat.tab.l[[c.fil]][cols$D2O]
  # Total AUCs for the D2O data (in duplicate) of the current dataset [cols$D2O-named num vector]
   aux.dat.norms <- dt.g*colSums(c.dat.df)
  # List normalized profiles (concatenated with index vector and time grid) [(H20,D20)-named list of df's]
   c.poll.df <- cbind(c.poll.df , sweep(c.dat.df,2,aux.dat.norms,'/')) # normalize and attach current D2O control data [df]
 }

# Extract and record in poll.prof.df the pollutant profile [(N.t)x(idx,x,y) df]
 aux.dat.cols <- which(!(colnames(c.poll.df) %in% c('idx','x'))) # data-holding cols in c.poll.df (integer(0) if fils.poll.v = NULL) [nat vector]
 poll.prof.df <- data.frame( # average normalized profiles (y-col = NaN if fils.poll.v = NULL)
   idx=c.poll.df$idx , x=c.poll.df$x , y=rowMeans(c.poll.df[,aux.dat.cols]) )
 if(all(is.nan(poll.prof.df$y))){poll.prof.df$y <- 0} # if no "pure pollutant" profiles (fils.poll.v = NULL), then make y-col zero
 poll.prof.df$y <- mollify( data.frame('x'=poll.prof.df$x , 'y'=poll.prof.df$y) )$y # mollify and renormalize pollutant profile
 poll.prof.df$y <- poll.prof.df$y - min(poll.prof.df$y) # remove offset

    print(sprintf('Tabulation of pollutant profile complete [%s]',organ))