####################################################################################################
# 'presets.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This routine prepares and stores the datasets corresponding to the organ set in settings.R (and
# stored in the global environment). Additionally, it defines or initializes a number of data
# structures feeding into our analysis.
#
# For the data structures, see index below (and inline for omissions). For the data, the module
# originally loads the time series from times.csv (in the data folder of the current organ) and then
# truncates it to the rough interval (cf. t.rough.lims.l in settings.R). The time series forms a
# regular grid over that interval but skips a few (approximately 1-2%)  of its values. We add these
# manually, obtaining a time grid t.g with step dt.g over [t.rough.lims$L , t.rough.lims$R].
#
# The profile data are initially stored in dat.tab.l, a list whose elements are df's. Each element
# is named after a dataset filename (i.e. an element of fils.v), has dat.cols as column names (see
# settings.R) and reports (in duplicate) the H2O and D2O profile counts over the original time
# series (stored in times.csv). Within this module, each such df is truncated to the rough interval
# for the current organ, as we don't use any information outside it. Missing profile values at the
# newly added time gridpoints ('dataset repair') are interpolated from the immediate neighbors. This
# way, dat.tab.l remains a dataset-named list, each elemenf of which is a df reporting (repaired)
# profile counts (columns arranged as in dat.cols) over the (regular) timegrid t.g.
#
#   INPUTS
# None.
#
#   OUTPUTS
# * dat.fld = path to folder containing the datasets for the current organ [string]
# * weeks.df = (raw) df holding the data in weeks.csv of the data folder [(Filename,ID,Weeks)-col df]
# * weeks.v = unique week numbers represented in weeks.df [num vector]
# * weeks.l = list of dataset filenames per week number [weeks.v-named list]
# * fils.v = (ordered) dataset filenames in data folder; the filenames include a '.csv' ending but
#            not the folder path - that is held by dat.fld [string vector]
# * plot.fils = df holding the (full) filepaths to the profile plots (plotted by plot.prof.R);
#               each row contains the 'data' column (name of .csv file) and the 'file' column
#               (filepath for that dataset) [length(fils.v)x(data,file) df]
# * t.g = uniform grid over the 'rough interval' [t.rough.lims$L,t.rough.lims$R] (cf. settings.R) [num vector]
# * dt.g = timestep for grid t.g [scalar]
# * N.t = cardinality of t.g [natural]
# * idx.g = sequence of consecutive naturals indexing the datasets (df's of dat.tab.l) [num vector]
# * dat.tab.l = list of datasets, named after the dataset files (i.e. fils.v), with each element
#                reporting the dataset (profile counts) of that file; the dataset is truncated to
#                the rough interval [t.rough.lims$L,t.rough.lims$R], i.e. all entries outside that
#                interval are discarded; missing values are interpolated linearly from their
#                immediate neighbors [(fils.v)-list of N.tx(idx,x,dat.cols) df's]
# t.rough.idxs = indices of the rough interval boundaries in the original (uncropped) t.g [(L,R) df]
#
#   DEPENDENCIES
# <- settings.R
# -> all other modules in folder 'Modules'
####################################################################################################
# PROCESS FOLDER AND FILE NAMES -->>

 # Path to data folder (containing the datasets for current organ) [string]
  dat.fld <- paste('../Results/',organ,'/CSV_data/',sep='')

 # Week numbers
  # Store weeks.csv in weeks.df [(Filename,ID,Weeks)-col df]
   weeks.df <- read.csv(paste(dat.fld,'weeks.csv',sep=''))
  # Identify and sort unique week numbers present in weeks.df [num vector]
   weeks.v <- sort(unique(weeks.df[!is.na(weeks.df$Weeks),'Weeks'])) # week numbers (ordered and w/o repetitions)
  # List of datasets (filenames) per (unique) week number [weeks.v-named list]
   weeks.l <- setNames(lapply(weeks.v , function(arg){as.character(weeks.df[!is.na(weeks.df$Weeks) & weeks.df$Weeks==arg,'Filename'])}) , weeks.v)

 # Order dataset filenames, then store them in fils.v [string vector]
  # List filenames of datasets present in the data folder
   fils.v <- list.files(dat.fld) # list all files in current dir... [string vector]
     fils.v <- fils.v[fils.v!='times.csv' & fils.v!='weeks.csv'] # ...then exclude times.csv and weeks.csv [subset of fils.v]
  # Reorder filenames of datasets by ascending numeral order
   # Make vector of dataset numerals in all data filenames [int vector]
    aux.length.organ <- length(unlist(strsplit(organ,''))) # no. of characters in the current organ name
    aux.fil.nos.v <- sapply(fils.v , function(arg){unlist(strsplit(arg,'[.]'))[[1]]}) # remove .csv ending from filenames [string vector]
    aux.fil.nos.v <- sapply(aux.fil.nos.v , function(arg){unlist(strsplit(arg,''))[-(1:(aux.length.organ+1))]}) # retain file numeral [string vector list]
    aux.fil.nos.v <- sapply(aux.fil.nos.v , function(arg){as.integer(paste(arg,collapse=''))}) # rejoin split numerals [int vector]
   # Make vector of indices (numerals) in all data filenames [int vector]
    fils.v <- fils.v[sort(aux.fil.nos.v , index.return=T)$ix] # reorder filename vector by ascending dataset numeral [permutation of fils.v]
   # Check consistency between fils.v and the filenames reported in weeks.csv
    if(any(fils.v != as.character(weeks.df$Filename))){stop(sprintf('Mismatch between files in %s folder and in weeks.csv$Filename [presets]',organ))}
    N.fils <- length(fils.v) # no. of datasets for current organ [nat]

 # Record control filenames (datasets) and their positions in fils.v for current organ [(fil,idx), (string,num) df]
  ctrl.fils.df <- data.frame(fil=fils.v[weeks.df$ID=='CTRL'] , idx=which(fils.v %in% fils.v[weeks.df$ID=='CTRL']) , stringsAsFactors=F)

 # Safety check: at least one control set is present and middle controls are absent (at most one index jump)
  aux.jmps <- which(diff(ctrl.fils.df$idx)>1) # index past which ctrl.fils.df$idx has jump > 1 [nat vector]
  if(!((1%in%ctrl.fils.df$idx) | (length(fils.v) %in% ctrl.fils.df$idx)) | length(aux.jmps)>1)
   {stop(sprintf('ABORT - check controls for %s [presets]',organ))}

 # Store start/end control sub-df's in a (start,end)-list; if one (start or end) set is empty, b/c the corresponding controls are absent,
 # then we copy into it the other, non-empty vector (resulting in identical start/end controls) [(start,end)-named list of ctrl.fils.df sub-df's]
  ctrl.fils.split.l <- list(start=c() , end=c()) # preallocate list holding start/end control datasets [(start,end)-named list]
  if(length(aux.jmps)==1)
  { # if a (single) jump exists, then both start and end control sets present
   ctrl.fils.split.l[['start']] <- head(ctrl.fils.df,aux.jmps)
   ctrl.fils.split.l[['end']] <- tail(ctrl.fils.df,-aux.jmps)
  }else
  { # if no jumps exist, then assign existing control set to both start/end control sets
   ctrl.fils.split.l[['start']] <- ctrl.fils.df
   ctrl.fils.split.l[['end']] <- ctrl.fils.df
  }

 # List non-pure-pollutant filenames (datasets) for current organ [string vector]
  dat.fils.df <- data.frame(fil=setdiff(fils.v,fils.poll.l[[organ]]) , idx=which(fils.v %in% setdiff(fils.v,fils.poll.l[[organ]])) , stringsAsFactors=F)

 # Make df holding plot pathways [length(fils.v)x(file,combi) df]
  plot.fils <- data.frame(plot=matrix(NA,nrow=length(fils.v)) , combi=matrix(NA,nrow=length(fils.v)) , row.names=fils.v) # preallocate df
  for(c.fil in fils.v)
  {
   plot.fils[c.fil,'plot'] <- paste('../Results/',organ,'/Plots/',strsplit(c.fil,'[.]')[[1]][1],'.pdf',sep='')
   plot.fils[c.fil,'combi'] <- paste('../Results/',organ,'/Combiplots/',strsplit(c.fil,'[.]')[[1]][1],'.pdf',sep='')
  }

# Organ-specific Width for noise-robust detection of peaks (see settings.R)
 N.w <- N.w.l[[organ]]

# PROCESS FOLDER AND FILE NAMES <<--
####################################################################################################
# INITIALIZE DATA STRUCTURES -->>

 # Make df holding start/end pure profiles [(start,end)-named list of dfs]
  pure.prof.l <- setNames(list(NA,NA) , names(ctrl.fils.split.l))

 # NA-preallocate list of df's holding left/right peak locations/boundaries and valley per start/end
 # pure profiles of current organ [(start,end)-list of (H2O,D2O)-lists (L,Lmax,M,Rmax,R)x(x,idx) NA df's]
  aux.v <- c(L=NA,Lmax=NA,M=NA,Rmax=NA,R=NA) # col template [(L,M,R) NA-vector]
  aux.df <- setNames(data.frame(aux.v , aux.v) , c('idx','x')) # individual list entry [(L,M,R)x(H2O,D2O) df]
  pure.loc.l <- list(start=list(H2O=aux.df,D2O=aux.df) , end=list(H2O=aux.df,D2O=aux.df))

 # NA-preallocate df's holding left/middle/right cutoffs of non-pollutant profiles [dat.fils.df-list of (idx,x)-lists of (L,M,R)x(dat.cols) df]
  prof.loc.l <- setNames(vector('list',length(dat.fils.df$fil)) , dat.fils.df$fil) # declare prof.loc.l [dat.fils.df$fil-named NULL list]
  loc.names <- c('L','M','R') # names of cutoffs [string vector]
  aux.df <- setNames(data.frame(matrix(NA,ncol=length(dat.cols),nrow=length(loc.names)),row.names=loc.names) , dat.cols) # [(loc.names)x(dat.cols)) df]
  for(c.fil in dat.fils.df$fil){prof.loc.l[[c.fil]] <- list(idx=aux.df , x=aux.df)} # preallocation

 # NA-preallocate df's holding position (start/end) of best matching pure profiles (per dataset), as well as L/M/R cutoffs for these profiles derived
 # by pushing fwd dataset cutoffs [dat.fils.df-list of (idx,x,orig)-lists of (L,M,R)x(dat.cols) df (idx-/x-elems) and 1x(dat.cols) df (orig-elem)]
  pure.crop.l <- setNames(vector('list',length(dat.fils.df$fil)) , dat.fils.df$fil)
  for(c.fil in dat.fils.df$fil){pure.crop.l[[c.fil]] <- list( idx=aux.df , x=aux.df ,
    orig=setNames(data.frame(matrix(NA,nrow=1,ncol=length(dat.cols))) , dat.cols)  )} # preallocation

 # NA-preallocate df's holding valley-derived cutoff candidates of non-pollutant profiles [dat.fils.df-list of (idx,x)-lists of (L,R)x(dat.cols) df]
  prof.valley.cand.l <- prof.loc.l
  for(c.fil in dat.fils.df$fil)
  { # remove M-row from df's
   prof.valley.cand.l[[c.fil]]$idx <- prof.valley.cand.l[[c.fil]]$idx[-which(rownames(prof.valley.cand.l[[c.fil]]$idx)=='M') , ]
   prof.valley.cand.l[[c.fil]]$x <- prof.valley.cand.l[[c.fil]]$x[-which(rownames(prof.valley.cand.l[[c.fil]]$x)=='M') , ]
  }

 # NA-preallocate df's holding pure-derived cutoff candidates of non-pollutant profiles [dat.fils.df-list of (idx,x)-lists of (L,R)x(dat.cols) df]
  prof.pure.cand.l <- prof.valley.cand.l

 # Preallocate list of df's holding non-pollutant dataset counts [fils.v-named NA list]
  cnt.l <- setNames(vector('list',length(dat.fils.df$fil)) , dat.fils.df$fil) # declare cnt.l [dat.fils.df$fil-named NULL list]
  cnt.l <- lapply(cnt.l,function(arg){arg <- NA}) # declare all list elements equal to NA [fils.v-named list of NA elements]

 # Preallocate list of df's holding repaired dataset counts [fils.v-named NA list]
  cnt.fix.l <- cnt.l

 # Preallocate df holding dataset skewness ratios (left interval counts over entire interval counts) [(dat.fils.df$fil)x(dat.cols) NA df]
  SR.df <- setNames(data.frame(matrix(NA,length(dat.fils.df$fil),length(dat.cols)) , row.names=dat.fils.df$fil) , dat.cols)

 # Preallocate list holding dataset ratios and skewness ratios (in duplicate, as well as their means) [(ratio,SR)-named empty list]
  ratio.prp.l <- list(ratio=c() , prp=c()) # declare ratio.prp.l [empty list]
  ratio.prp.l$prp <- data.frame(week=c() , A=c() , B=c() , AB=c()) # prp-element of ratio.prp.l

 # Preallocate list holding (original and repaired) week-averaged (aggregate) ratios and signal proportions [(orig,fix)-list of (week,L,R,T,prp)-col df]
  ratio.SR.week.aggr.l <- list(
    orig = data.frame(week=weeks.v,L=NA,R=NA,T=NA,prp=NA) , # df for original ratios
    fixd = data.frame(week=weeks.v,L=NA,R=NA,T=NA,prp=NA) ) # df for repaired ratios

 # Initialize df holding prp-values, i.e. proportions of signal (per dataset; in duplicate) [(dat.fils.df$fil)x(cols$D2O) df]
  c.prp.df <- setNames(data.frame(matrix(NA,nrow=length(dat.fils.df$fil),ncol=2) , row.names=dat.fils.df$fil) , cols$D2O)

# INITIALIZE DATA STRUCTURES <<--
####################################################################################################
# READ AND STORE TIME SERIES AND PROFILE COUNTS -->>

 # Store time grid in t.g [num vector]
  t.g <- read.csv(paste(dat.fld,'times.csv',sep=''))[,1] # timepoints [num vector]
  N.t <- length(t.g) # initial (i.e. pre-cropping to rough interval) no. of timepoints [natural]

 # Make df holding original datasets (count profiles) for current organ [(fils.v)-list of length(t.g)x(dat.cols) df's]
  dat.tab.l <- setNames(vector('list',length(fils.v)) , fils.v) # declare dat.tab.l  [fils.v-named NULL list]
  for(c.fil in fils.v)
  { # iterate over datasets in data folder
   # Read dataset, remove NA rows [(N.t)x(dat.cols) df]
    dat.tab.l[[c.fil]] <- read.csv(paste(dat.fld,c.fil,sep='')) # read dataset
    dat.tab.l[[c.fil]] <- dat.tab.l[[c.fil]][apply(dat.tab.l[[c.fil]],1,function(arg){all(!is.na(arg))}),] # crop NA rows
   # Rename columns (cf. settings.R) after a sanity check
    if(dim(dat.tab.l[[c.fil]])[1]!=N.t){stop(sprintf('%s-data length differs from time series length',c.fil))}else{names(dat.tab.l[[c.fil]]) <- dat.cols}
  } # end of iteration over datasets in data folder

# READ AND STORE TIME SERIES AND PROFILE COUNTS <<--
####################################################################################################
# TRUNCATE AND REPAIR TIME SERIES AND PROFILE COUNTS -->>

 # Left/right limits of rough interval (in time domain) for current organ [(L,R)-col df]
  t.rough.lims <- t.rough.lims.l[[organ]]

 # Truncate time series to rough interval, then add missing gridpoints [num vector]
  t.rough.idxs <- data.frame('L'=min(which(t.rough.lims$L <= t.g)) , 'R'=max(which(t.g <= t.rough.lims$R))) # indices of t.rough.lims in t.g [(L,R) df]
  t.g <- t.g[t.rough.idxs$L:t.rough.idxs$R] # truncate t.g to rough interval
  dt.g <- min(diff(t.g)) # timestep [pos scalar]
  t.g.miss <- 1+which(diff(t.g , lag=1 , differences=2) > 1E3*.Machine$double.eps) # missing timepoints (indices) [nat vector]
  t.g.adds <- t.g[t.g.miss]/2 + t.g[t.g.miss+1]/2 # missing timepoints (values) [num vector]
  for(c.n in length(t.g.miss):1){t.g <- append(x=t.g , value=t.g.adds[c.n] , after=t.g.miss[c.n])} # insert missing timepoints [num vector]
  N.t <- length(t.g) # update no. of timepoints [natural]
    if(any(diff(t.g , lag=1 , differences=2) > 1E3*.Machine$double.eps)){stop('Missing timepoint problem persists [presets]')} # sanity check

 # Truncate data to rough interval [t.rough.lims$L,t.rough.lims$R], then repair them [(fils.v)-list of length(t.g)x(idx,x,dat.cols) df's]
  for(c.fil in fils.v)
  { # iterate over datasets in data folder
   # Truncate dataset
    dat.tab.l[[c.fil]] <- dat.tab.l[[c.fil]][t.rough.idxs$L:t.rough.idxs$R,]
   # Interpolate and append values at missing timepoints [N.t-long df]
    c.dat.df.adds <- round(dat.tab.l[[c.fil]][t.g.miss,]/2 + dat.tab.l[[c.fil]][t.g.miss+1,]/2) # interpolated missing counts [num vector]
    for(c.n in length(t.g.miss):1)
    { # insert interpolated missing counts (end-to-start, to retain crucial indexing during iteration)
     dat.tab.l[[c.fil]] <- rbind( head(dat.tab.l[[c.fil]],t.g.miss[c.n]) , c.dat.df.adds[c.n,] , tail(dat.tab.l[[c.fil]],-t.g.miss[c.n]) )
    }
   # Subtract constant offset (minimum per column) [(N.t)x(dat.cols) df]
    dat.tab.l[[c.fil]] <- sweep(dat.tab.l[[c.fil]] , 2 , sapply(dat.tab.l[[c.fil]],min) , '-') # remove background
   # Reindex df
    rownames(dat.tab.l[[c.fil]]) <- seq(from=as.integer(rownames(dat.tab.l[[c.fil]])[1]),length.out=N.t,by=1)
      if(dim(dat.tab.l[[c.fil]])[1]!=N.t){stop(sprintf('Length of %s data differs from length of time series',c.fil))} # sanity check
   # Add time (x) and index (idx) column, then bring them forth (helps visually)
    dat.tab.l[[c.fil]]$x <- t.g
    dat.tab.l[[c.fil]]$idx <- as.integer(rownames(dat.tab.l[[c.fil]]))
    dat.tab.l[[c.fil]] <- dat.tab.l[[c.fil]][,c('idx','x',names(dat.tab.l[[c.fil]])[!(names(dat.tab.l[[c.fil]]) %in% c('idx','x'))])] # rearrange cols
  } # end of iteration over datasets in data folder

 # Check all datasets have identical index and time cols, then define global index vector
  for(c.fil in fils.v){if(any(dat.tab.l[[c.fil]]$idx!=dat.tab.l[[fils.v[1]]]$idx)){stop('ABORT: check file indexing [presets.R]')}}
  for(c.fil in fils.v){if(any(dat.tab.l[[c.fil]]$x!=dat.tab.l[[fils.v[1]]]$x)){stop('ABORT: check file timegriding [presets.R]')}}
  idx.g <- dat.tab.l[[fils.v[1]]]$idx # sequence of consecutive natural numbers indexing the datasets [N.t-long int vector]

# TRUNCATE AND REPAIR TIME SERIES AND PROFILE COUNTS <<--
####################################################################################################

    print(sprintf('Presets complete [%s]',organ))