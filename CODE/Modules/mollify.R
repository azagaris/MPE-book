####################################################################################################
# 'mollify.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This routine defines the function mollify, which uses a Gaussian kernel to mollify profiles passed
# to it (as data frames with columns 'x', 'y' and possibly others). The function is largely copied
# from the earlier function mollifier.R, which I wrote for CHROME (Ad's chromatography project).
# Since mollification perturbs the AUC slightly due to boundary effects, profiles are renormalized
# (post-processed) to recover their original AUC (equivalently, the total sum of their y-entries).
# All columns other than the y-col are left untouched. In particular, the output df has the same
# structure as the input df, with columns others than y returned in their original state.
#
#   INPUTS
# * prof.in = a df with x- and y-cols and possibly (but not necesssarily) others; the x-column
#             should form a regular grid [(x,y,...)-col df]
#
#   OUTPUTS
# * prof.mol = a df with the same (type of) columns as prof.in; all cols other than y are identical
#              to the corresponding ones of prof.in, whereas the y-col holds the mollifed y-values
#                                                                                 [(x,y,...)-col df]
#
#   DEPENDENCIES
# <- settings.R
# -> extrema.locs.R [if smoothen=TRUE in extrema.locs.R]
# -> pure.tab.R
# -> poll.tab.R
####################################################################################################
mollify <- function(prof.in){

 # Check input is (x,y)-col df with regular x-grid
  if(!is.data.frame(prof.in) || !all(c('x','y') %in% names(prof.in))){stop('ABORT: check input type [mollify]')} # abort if input is of wrong type
  if(any(diff(diff(prof.in$x,1),1) > peak.thr*max(diff(prof.in$x,1)))){stop('ABORT: x-grid irregular [mollify]')} # abort if x-grid is (too) irregular

 # Dataset length and AUC [scalars]
  N.t.in <- dim(prof.in)[1]
  AUC.in <- dt.g*sum(prof.in$y)

 # Tabulate mollifier coefficients [(2L+1)-dim symmetric vector]
  mollify.coeffs <- exp(-seq(from=-L,to=L,by=1)^2/(2*sgm^2)) # (unnormalized) Gaussian coefficients
  mollify.coeffs <- mollify.coeffs/sum(mollify.coeffs) # normalize coefficients to unit sum

 # Extend (horizontally) tabulated profile prof.in beyond endpoints
  prof.ext <- prof.in[,c('x','y')] # initialize extended tabulation [(x,y)-named df]
  prof.ext.L <- data.frame( # left extension
    'x' = seq(from = prof.in$x[1]-L*dt.g , by = dt.g , length.out = L) , # left-extend grid 
    'y' = rep(prof.in$y[1] , L) ) # left-extend tabulated function horizontally
  prof.ext.R <- data.frame( # right extension
    'x' = seq(from = prof.in$x[N.t.in]+dt.g , by = dt.g , length.out = L) , # right-extend grid
    'y' = rep(prof.in$y[N.t.in],L) ) # right-extend tabulated function horizontally
  prof.ext <- rbind(prof.ext.L,prof.ext,prof.ext.R) # append left and right extensions to tabulation

 # Mollify profile and store in prof.mol [(N.t.in)x(x,y) df]
  prof.mol <- prof.in # preallocate output with input (retains possible extra columns)
  C.mlf <- matrix(0,N.t.in+2*L,N.t.in+2*L) # preallocate mollification matrix
  for(idx in 1:N.t.in){ C.mlf[idx+L,idx:(idx+2*L)] <- mollify.coeffs } # implant coefficients in idx-th row
  prof.mol$y <- ( C.mlf %*% as.vector(t(prof.ext$y)) )[(L+1):(N.t.in+L)] # convolution

 # Renormalization to original AUC (unless AUC=0) [(N.t.in)x(x,y) df]
  if(AUC.in>0){prof.mol$y <- prof.mol$y*(AUC.in)/sum(dt.g*prof.mol$y)}

    return(prof.mol)
} # end of function