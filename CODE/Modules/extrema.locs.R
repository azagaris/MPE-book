####################################################################################################
# 'extrema.locs.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This function identifies (essential) maxima and minima of an input profile prof.in. 'Essential'
# means that they are likely to be extrema of the underlying smooth profile, rather than be due to
# noise on top of that profile. To assist with discarding spurious (i.e. non-essential) extrema, the
# function includes the option of ameliorating noise by mollification of the input profile. This
# option should be used with moderation, as mollification adversely affects the function near the
# boundary (where we extend it as constant).
#
# The input profile contains an index column ('idx'), timegrid column ('x') and count column ('y').
# Identification proceeds by finding all local maxima and eliminating those due to noise. In more
# detail, we first find all times with counts no smaller than those of the point's N.w-many left and
# N.w-many right neighbors (N.w is set in settings.R). An easy calculation establishes that, under
# the assumption that the profile is constant and the noise symmetric (i.e. a point is above/below
# an immediate neighbor with equal probabilities), pure noise will exhibit no peaks (as defined
# above) with probability (1 - 1/4^N.w)^(N.t-2*N.w) - here, we exclude N.w-wide buffers. Since
# N.w << N.t (buffer << time domain) for the profiles we consider, we find N.w ~ log_4(eps/N.t) for
# probability 1-eps of finding no peaks in such a signal. For eps ~1/100, we thus obtain N.w ~ 9.
# Identification of the minima proceeds in much the same way.
#
# Side notes: the module discards peaks peaking lower than (100*peak.thr)% of the global maximum
# (highest peak), as they are deemed unsubstantial (i.e. due to noise); note that peak.thr is a
# threshold set in settings.R. This practice is controversial, as the presence of contaminant in
# highly polluted profiles can elevate hopelessly the global maximum, which in turn leadd to weeding
# out substantial peaks. Similarly, only interpeak minima (between the lowest and highest identified
# peaks) were originally retained and the rest discarded. This practice was abandoned, since the
# routine is now (also) used to set cutoffs for dataset profiles, and minima below/above the
# secondary/principal peak are relevant there. Lastly, a final sanity check that maxima and minima
# were intertwined was also removed, as it offered nothing and scenarios where successive minima
# occur without a maximum in between can actually happen.
#
#   INPUTS
# * prof.in = a df with x- and y-cols and possibly (but not necesssarily) others [(x,y,...)-col df]
#
#   OUTPUTS
# * out.l = a list with max and min elements; each such element is a list itself, with loc and idx
#           elements; the loc element lists the x-locations of the min or max extreme, depending on
#           parent element; the idx element lists the indices corresponding to these extrema, as
#           these are listed in the input's idx-column; if such a column is missing, then the idx
#           element is NA [(max,min)-list of (idx,x)-df's]
#
#   DEPENDENCIES
# <- settings.R
# <- mollify.R [if smoothen=TRUE]
# -> pure.tab.R
# -> poll.tab.R
####################################################################################################
extrema.locs <- function(prof.in,smoothen=FALSE){

 # Check input is (x,y)-col df
  if(!is.data.frame(prof.in) || !all(c('x','y') %in% names(prof.in))){stop('ABORT: check input type [extrema.locs]')}
 # Dataset length [natural]
  N.t.in <- dim(prof.in)[1]
 # Set switch indicating whether an index-column is present [logi]
  idx.swtc <- 'idx' %in% names(prof.in)
 # Define index column, if absent
  if(!idx.swtc){prof.in$idx <- 1:N.t.in}
 # Mollify input profile, if smoothening switch turned on
  if(smoothen){prof.in <- mollify(prof.in)}
 # Preallocate output [(max,min)-list]
  out.l <- list(max=NA , min=NA)

 # Identify peaks (maxima)
  # Store into aux.df the idx-col, y-col (excluding N.w-wide boundary buffers) and progressive left/right y-col shifts [(N.t.in-2*N.w)x(2+2*N.w) df]
   aux.df <- data.frame(matrix( NA , nrow=N.t.in-2*N.w , ncol=2+2*N.w , dimnames=list(c(),c('idx','y',1:(2*N.w))) )) # preallocate aux.df
   aux.df[,'idx'] <- prof.in$idx[(1+N.w):(N.t.in-N.w)] # 1st col: row no's in prof.in of data (excluding boundary buffers)
   aux.df[,'y'] <- prof.in$y[(1+N.w):(N.t.in-N.w)] # 2nd col: original data (i.e. y-col excluding boundary buffers)
   for(c.n in 1:N.w)
   { # iterate over neighbor orders (i.e. degree of closeness)
    aux.df[,2+c.n] <- prof.in$y[(1+N.w-c.n):(N.t.in-N.w-c.n)] # (c.n+2)-nd col: c.n-th closest left neighbors
    aux.df[,2+c.n+N.w] <- prof.in$y[(1+N.w+c.n):(N.t.in-N.w+c.n)] # (N.w+c.n+2)-nd col: c.n-th closest right neighbors
   }
  # Identify indices of points no smaller than their N.w-many left/right closest neighbors
   aux.idxs <- aux.df[rowSums(aux.df[,'y'] < aux.df[,!(names(aux.df) %in% c('idx','y'))])==0 , 'idx'] # indices of local maxima [subset of aux.df$idx]
  # Weed out peaks (maxima) being too small (i.e. below peak.thr-times the global maximum)
   aux.idxs <- aux.idxs[ which(aux.df[aux.df$idx %in% aux.idxs,'y'] > peak.thr*max(aux.df$y)) ] # only retain large enough peaks
   aux.idxs <- aux.idxs[ sort(aux.idxs,decreasing=FALSE,index.return=T)$ix ] # order indices (ascending)
   aux.idxs <- aux.idxs[ which(c(diff(aux.idxs,1),Inf)!=1) ] # break ties (retain one of two adjacent peaks) [vector]
  # Record in output list the locations (times) and indices of substantial peaks
   out.l$max <- list(idx=aux.idxs , x=prof.in[prof.in$idx %in% aux.idxs,'x']) # store locations and indices [(x,idx)-list]
   if(!idx.swtc){out.l$max$idx <- NA} # if indes col absent from input df, do not record indices in output
   if(length(out.l$max$x)<2){stop('ABORT: less than two peaks found [extrema.locs]')} # verify at least two peaks found

 # Identify valleys (minima)
  # Identify indices of points no larger than their N.w-many left/right closest neighbors
   aux.idxs <- aux.df[rowSums(aux.df[,'y'] > aux.df[,!(names(aux.df) %in% c('idx','y'))])==0 , 'idx'] # indices of local minima [subset of aux.df$idx]
   aux.idxs <- aux.idxs[ sort(aux.idxs,decreasing=FALSE,index.return=T)$ix ] # order indices (ascending)
   if(length(aux.idxs)>1){aux.idxs <- aux.idxs[ which(c(diff(aux.idxs,1),Inf)!=1) ]} # break ties (retain one of two adjacent peaks) [vector]
  # Record in output list the locations (times) and indices of substantial peaks
   out.l$min <- list(idx=aux.idxs , x=prof.in[prof.in$idx %in% aux.idxs,'x']) # store locations and indices [(x,idx)-list]
   if(!idx.swtc){out.l$min$idx <- NA} # if indes col absent from input df, do not record indices in output

 # Turn (max,min)-list of (x,idx)-lists into (max,min)-list of (x,idx)-df's
  for(c.extr in c('max','min')){out.l[[c.extr]] <- data.frame(out.l[[c.extr]])}

     return(out.l)

 } # end of function