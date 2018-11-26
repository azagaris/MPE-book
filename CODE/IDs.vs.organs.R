####################################################################################################
# 'IDs.vs.organs.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This module inventories available datasets, first assigning them to particular IDs (i.e. goats)
# and organs (i.e. cell types) and then coloring them depending on whether they are pure pollutant
# profiles or not.

# Specifically, the current module collects all non-control datasets from non-blood organs (meaning
# the dataset ensembles Blood_CD11b_14_Batch_1 and Blood_CD21 are entirely excluded, as are datasets
# marked CTRL in the remaining organs). Each dataset has (in weeks.csv, within the organ-specific
# folder under the folder 'Results') a full ID of the form 'XXX T=Y'. Here, XXX is a 3-digit numeral
# (goat ID) and Y an integer. For all organs other than Ileum_CD14, Y equals the entry in the 'Weeks'
# column of the same .csv file; for Ileum_CD14, Y equals that entry plus two. (A few deviations from
# this rule were identified using the module and passed on to Lars, as they are probably mistakes
# affecting the ever0important 'Weeks' entry; by the time you, dear reader, are reading this, they
# might be present or not depending on whether they will have been corrected. Similarly, I suggested
# that Y be uniformly equal to the 'Weeks' entry for all non-blood organs, if that makes sense, in
# which case it could also be abolished; this suggestion might also have been implemented or not.)
#
# Blood samples cannot be subjected to that analysis, since blood was repeatedly drawn from specific
# goats. (Actual) organ samples, though, demanded that goats be sacificed, meaning that each such
# dataset only corresponds to a single 'Weeks' entry and value of Y.
#
# Once all full IDs have been inventoried, the 'T=Y' part is removed to yield numerical (i.e. goat)
# IDs of the form 'XXX'. Because of different conventions (see above) on Y-values and Weeks-entries
# for Ileum_CD14 and the rest of (non-blood) organs, there are roughly twice as many full IDs as
# numerical IDs (by the time I'm writing this anyway). These numerical IDs and the organs are paired
# in the df ID.organ.df, whose entries are character strings. The entry 'gray' at the (i,j) position
# signifies that the i-th goat did not yield a sample of the j-th organ. The entry 'red' at that
# position signifies, instead, that the i-th goat did yield a sample of the j-th organ, and that
# profile was chosen (in settings.R) as a pure pollutant profile. Finally, the entry 'blue' at that
# position signified availability of a j-th organ sample from the i-th goat which, additionally, is
# not a pure pollutant profile.
#
#
#   INPUTS
# None.
#
#   OUTPUTS
# * nb.organs.poll = nb.organs-subset listing organs w/ pure-pollutant profiles [string vector]
#
# * nb.organs.no.poll = nb.organs-subset listing organs w/o pure-pollutant profiles [string vector]
#
# * all.IDs.v = all full IDs (of the form 'XXX T=Y') [string vector]
#
# * all.IDs.l = all full IDs split per organ [nb.organs-named list of string vectors]
#
# * all.IDs.num.v = all numerical IDs (of the form 'XXX') [string vector]
#
# * all.IDs.num.l = all numerical IDs split per organ [nb.organs-named list of string vectors]
#
# * poll.IDs.v = all full IDs yielding (for some organ) pure-pollutants [string vector]
#
# * poll.IDs.l = all full IDs yielding (for some organ) pure-pollutants; split per organ
#                                                      [nb.organs.poll-named list of string vectors]
#
# * poll.IDs.num.v = all numerical IDs yielding (for some organ) pure-pollutants [string vector]
#
# * poll.IDs.num.l = all numerical IDs yielding (for some organ) pure-pollutants; split per organ
#                                                           [nb.organs-named list of string vectors]
#
# * ID.organ.df = a df with entries either 'gray', 'red' or 'blue', chosen as described above
#                 depending on the availability and quality of specific organ samples from specific
#.                goats [(all.IDs.v)x(nb.organs) df]
#
#   DEPENDENCIES
# <- settings.R
####################################################################################################
source('Modules/settings.R')

# Split organs into organs with/without pure-pollutant profiles [string vectors]
 nb.organs.poll <- intersect(nb.organs,names(fils.poll.l)) # ...with...
 nb.organs.no.poll <- setdiff(nb.organs,names(fils.poll.l)) # ...without...

# Preallocate lists to hold, per organ, all IDs (all.IDs.l) and pure-pollutant IDs (poll.IDs.l) [lists]
 all.IDs.l <- list() # declare all.IDs.l [empty list]
   length(all.IDs.l) <- length(nb.organs) # initialize all elements of all.IDs.l by NULL [nb.organs-named NULL list]
   names(all.IDs.l) <- nb.organs # name list elements by non-blood organs
 poll.IDs.l <- list() # declare all.IDs.l [empty list]
   length(poll.IDs.l) <- length(nb.organs.poll) # initialize all elements of poll.IDs.l by NULL [nb.organs.poll-named NULL list]
   names(poll.IDs.l) <- nb.organs.poll # name list elements by non-blood organs having pure-pollutant profiles

# Record, per organ, all/pure-pollutant IDs [lists]
 for(c.organ in nb.organs)
 { # iterate over organs having pure-pollutant profiles
  aux.df <- read.csv(paste(paste('../Results/',c.organ,'/CSV_data/',sep='') , 'weeks.csv' , sep='')) # read current organ summary
  all.IDs.l[[c.organ]] <- as.character(aux.df$ID)[as.character(aux.df$ID) != 'CTRL'] # record non-control IDs [organs-named list]
  if(c.organ %in% nb.organs.poll) # record pure-pollutant IDs [organs.poll-named list]
   {poll.IDs.l[[c.organ]] <- as.character(aux.df[aux.df$Filename %in% fils.poll.l[[c.organ]] , 'ID'])}
 }

# Vectorize all/pure-pollutant IDs (i.e. take union over all non-blood organs) [string vectors]
 all.IDs.v <- sort(unique(unlist(all.IDs.l))) # all (unique) IDs including 'T=...'
 poll.IDs.v <- sort(unique(unlist(poll.IDs.l))) # pure-pollutant (unique) IDs including 'T=...'

# Remove ' T=...' from all/pure-pollutant ID lists and vectors (i.e. only retain the numeral in each entry)
 all.IDs.num.l <- lapply(all.IDs.l , function(arg){sapply(arg,function(arg.aux){strsplit(arg.aux,' ')[[1]][1]})})
 poll.IDs.num.l <- lapply(poll.IDs.l , function(arg){sapply(arg,function(arg.aux){strsplit(arg.aux,' ')[[1]][1]})})
 all.IDs.num.v <- sort(unique(sapply(all.IDs.v , function(arg){strsplit(arg,' ')[[1]][1]}))) # all (unique) IDs excluding 'T=...'
 poll.IDs.num.v <- sort(unique(sapply(poll.IDs.v , function(arg){strsplit(arg,' ')[[1]][1]}))) # pure-pollutant (unique) IDs excluding 'T=...'

# Check numerals within each non-blood organ are unique (each goat is sacrificed once and yields at most one sample per organ)
 if(!all(sapply(all.IDs.num.l , function(arg){ifelse(length(unique(arg))==length(arg),T,F)}))){stop('ABORT - duplicate num IDs in non-blood organs')} 

# Check, for each numerical ID, how many weeks does it correspond to
 organ.ID.week.l <- list() # declare organ.ID.week.l [empty list]
   length(organ.ID.week.l) <- length(all.IDs.num.v) # initialize all elements of organ.ID.week.l by NULL [all.IDs.num.v-named NULL list]
   names(organ.ID.week.l) <- all.IDs.num.v # name list elements by numerical IDs
 for(c.ID.num in all.IDs.num.v)
 { # iterate over numerical IDs
  aux.organs.l <- nb.organs[sapply(all.IDs.num.l , function(arg){c.ID.num %in% arg})] # organs exhibiting current ID [nb.organs-named list]
  for(c.organ in aux.organs.l)
  { # iterate over organs containing current numerical ID
   aux.df <- read.csv(paste(paste('../Results/',c.organ,'/CSV_data/',sep='') , 'weeks.csv' , sep='')) # read current organ summary
   c.ID <- names(all.IDs.num.l[[c.organ]][all.IDs.num.l[[c.organ]]==c.ID.num]) # unique full name assumed by current numerical ID within current organ
   c.week <- aux.df[aux.df$ID==c.ID,'Weeks'] # week of current ID within current organ
   organ.ID.week.l[[c.ID.num]] <- rbind(organ.ID.week.l[[c.ID.num]] , data.frame(num.ID=c.ID.num , organ=c.organ , full.ID=c.ID , week=c.week))
  } # end of iteration over organs containing current numerical ID
 } # end of iteration over numerical IDs

# ID-vs-organ-status df: black if ID absent from organ, blue/red if present and clean/pollutant [(all.IDs.num.v)x(nb.organs) df]
 ID.organ.df <- data.frame(matrix(NA , nrow=length(all.IDs.num.v) , ncol=length(nb.organs) , dimnames=list(all.IDs.num.v,nb.organs))) # preallocate df
 ID.organ.df[,] <- 'gray' # initialize df w/ gray
 for(c.organ in nb.organs)
  { # iterate over non-blood organs
   ID.organ.df[all.IDs.num.l[[c.organ]],c.organ] <- 'blue' # all present IDs gray -> blue
   if(c.organ %in% nb.organs.poll){ID.organ.df[poll.IDs.num.l[[c.organ]],c.organ] <- 'red'} # all pollutant IDs gray/blue -> red
  }

# Plot results
 plot(c()  ,c() , typ='p' , xlim=c(1/dim(ID.organ.df)[1],1) , ylim=c(1/dim(ID.organ.df)[2],1) , xlab='' , ylab='' , xaxt='n' , yaxt='n')
 for(c.j in 1:dim(ID.organ.df)[2]){lines(c(-1,1) , c(c.j,c.j)/dim(ID.organ.df)[2] , typ='l' , col='MintCream')} # plot vertical lines
 for(c.i in 1:dim(ID.organ.df)[1]){lines(c(c.i,c.i)/dim(ID.organ.df)[1] , c(-1,1) , typ='l' , col='MintCream')} # plot horizontal lines
 for(c.i in 1:dim(ID.organ.df)[1]) # plot points
 { # iterate over numerical IDs (x-values)
  for(c.j in 1:dim(ID.organ.df)[2])
  { # iterate over non-blood organs (y-values)
   lines(c.i/dim(ID.organ.df)[1] , (1+dim(ID.organ.df)[2]-c.j)/dim(ID.organ.df)[2] , typ='p' , pch=15 , cex=1.75 , col=ID.organ.df[c.i,c.j])
  }
 }
 axis(side=1, at=(1:dim(ID.organ.df)[1])/dim(ID.organ.df)[1] , labels=all.IDs.num.v , las=2 , cex.axis=0.75) # horizontal ticks (numerical IDs)
 axis(side=2, at=(dim(ID.organ.df)[2]:1)/dim(ID.organ.df)[2] , labels=nb.organs , las=2 , cex.axis=0.33) # vertical ticks (non-blood organs)