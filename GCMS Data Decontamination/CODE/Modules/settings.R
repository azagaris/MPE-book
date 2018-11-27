####################################################################################################
# 'settings.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This routine contains global settings used by the subroutines called by the wrapper main.R. The
# settings are split in tunable settings (e.g. organ name whose datasets will undergo analysis),
# which end-users can modify, and fixed settings; these last are reserved for programmers and should
# not modified by end-users.
#
#   INPUTS
# None.
#
#   OUTPUTS
# None.
#
#   DEPENDENCIES
# -> all modules in folder 'Modules'
####################################################################################################
# END-USER TUNABLE SETTINGS -->>

 # Name of organ to analyze in current run [organs-element]
  organ <- 'Jejunum_CD14'

 # Plot switch, controlling whether plots should be redrawn and written, in the Plots and Combiplots folders, or not
  plot.swtc <- TRUE#FALSE# # set plot switch [logi]

 # Write switches, controlling whether analysis outcomes should be rewritten, in the CSV_Results folder, or not
  write.swtc <- FALSE#TRUE# # set plot switch [logi]

 # Threshold for detection of substantial peaks: peaks < (peak.thr)*(global max) are discarded [scalar in (0,1)]
  peak.thr <- 5/100

 # Standard deviation of Gaussian mollifier (in timesteps), controlling kernel broadness
  sgm <- 3L # set sd [natural, e.g. 7L]

 # Half-length of mollifier support (in timesteps), controlling kernel cutoff
  L <- 3L*sgm # set half-length [natural multiple of sgm, e.g. 3L*sgm]

# END-USER TUNABLE SETTINGS <<--
####################################################################################################
# FIXED SETTINGS - DO NOT TOUCH UNLESS RESTRUCTURING CODEBASE -->>

 # Names of all organs to be analyzed [string vector]
  organs <- list.files(paste('../Results/',sep='')) # all organs for which there are datasets
  nb.organs <- setdiff(organs,c('Blood_CD11b_14_Batch_1','Blood_CD21')) # non-blood organs
    if(!(organ %in% organs)){stop('ABORT: check organ name [settings.R]')}

 # Names of various categories used in the scripts [string vectors]
  mols <- c('H2O','D2O') # labels differentiating undeuterated (H2O) and deuterated (D2O) profiles
  dupls <- c('A','B') # labels differentiating duplicates per sample
  poss <- c('start','end') # labels differentiating control sets (and thus also pure profiles)

 # Column names to apply to each dataset (mols-and-dupls combinations; order dictated by .cvs files) [string vector]
  dat.cols <- c('H2O.A','D2O.A','H2O.B','D2O.B') # datasets column names (renamed b/c not unique, when loaded from .csv files)

 # Split dat.cols into undeuterated (H2O)/deuterated (D2O) sets and duplicates (A/B) [(H2O,D2O,A,B)-named list of 2-elem string vectors]
  cols <- list(NA,NA,NA,NA) # preallocate columns of undeuterated (H2O) and deuterated (D2O) data [(H2O,D2O)-list]
    names(cols) <- c(mols,dupls)
  for(c.mol in mols){ cols[[c.mol]] <- dat.cols[sapply( dat.cols , function(arg){strsplit(arg,'[.]')[[1]][1]} ) == c.mol] }
  for(c.dupl in dupls){ cols[[c.dupl]] <- dat.cols[sapply( dat.cols , function(arg){strsplit(arg,'[.]')[[1]][2]} ) == c.dupl] }

 # Column names for the df storing counts (dat.cols-cols) and deuteration ratios per duplicate (dupls-cols) and averaged (last col) [string vector]
  ratio.cols <- c(dat.cols,dupls,paste(dupls,collapse='')) # string vector (H2O.A,D2O.A,H2O.B,D2O.B,A,B,AB)

 # Row names for the df's stored within each element of cnts.l [string vector]
  ratio.rows <- c('L','R','T')

 # Width for noise-robust detection of peaks: a data point counts as peak iff >= its N.w-many left/right neighbors [natural]
  N.w.l <- list( # see preamble of pure.prof.lims.R for an educated guess; finetuning aids weeding out spurious minima in dat.vs.pure.R
   Blood_CD11b_14_Batch_1 = 10L ,
   Blood_CD21 = 10L ,
   Ileum_CD14 = 15L ,
   Ileum_CD21 = 10L ,
   Jejunum_CD11b = 10L ,
   Jejunum_CD14 = 10L ,
   Jejunum_CD21 = 10L ,
   Lung_CD11b_14 = 10L ,
   Lung_Wash_CD21 = 10L ,
   Mes_LN_CD14 = 10L ,
   Mes_LN_CD21 = 10L ,
   Peritoneal_Wash_CD21 = 15L ,
   Prescap_LN_CD14 = 10L ,
   Prescap_LN_CD21 = 15L )
     if(!(all(names(N.w.l) %in% organs) & all(organs %in% names(N.w.l)))){stop('ABORT: check N.w.l names [settings.R]')}

 # Left/right limits of rough interval (in time domain) enclosing both peaks across all datasets [(L,R)-col df]
  t.rough.lims.l <- list(
   Blood_CD11b_14_Batch_1 = data.frame(L=5.05 , R=5.27) ,
   Blood_CD21 = data.frame(L=5.00 , R=5.20) ,
   Ileum_CD14 = data.frame(L=5.04 , R=5.25) ,
   Ileum_CD21 = data.frame(L=5.01 , R=5.20) ,
   Jejunum_CD11b = data.frame(L=5.03 , R=5.20) ,
   Jejunum_CD14 = data.frame(L=5.05 , R=5.25) ,
   Jejunum_CD21 = data.frame(L=5.00 , R=5.25) ,
   Lung_CD11b_14 = data.frame(L=5.00 , R=5.25) ,
   Lung_Wash_CD21 = data.frame(L=5.00 , R=5.20) ,
   Mes_LN_CD14 = data.frame(L=5.00 , R=5.20) ,
   Mes_LN_CD21 = data.frame(L=5.00 , R=5.25) ,
   Peritoneal_Wash_CD21 = data.frame(L=5.00 , R=5.20) ,
   Prescap_LN_CD14 = data.frame(L=5.00 , R=5.25) ,
   Prescap_LN_CD21 = data.frame(L=5.00 , R=5.25) )
     if(!(all(names(t.rough.lims.l) %in% organs) & all(organs %in% names(t.rough.lims.l)))){stop('ABORT: check t.rough.lims.l names [settings.R]')}

 # List holding indices (.csv file numbers) of pure-contaminant profiles per organ [organs-named list of string vectors]
  aux.l <- list(
   Blood_CD11b_14_Batch_1 = c(22:30) ,
   Blood_CD21 = NA ,
   Ileum_CD14 = c(16,38:39) ,
   Ileum_CD21 = c(7,10:12,14:19,21:22) ,
   Jejunum_CD11b = c(6) ,
   Jejunum_CD14 = c(15:17,24,34:36) ,
   Jejunum_CD21 = NA ,
   Lung_CD11b_14 = NA ,
   Lung_Wash_CD21 = c(9,12:19) ,
   Mes_LN_CD14 = NA ,
   Mes_LN_CD21 = NA ,
   Peritoneal_Wash_CD21 = c(8,15:17) ,
   Prescap_LN_CD14 = NA ,
   Prescap_LN_CD21 = c(22,28,30:35,39:46) )
     if(!(all(names(aux.l) %in% organs) & all(organs %in% names(aux.l)))){stop('ABORT: check aux.l names [settings.R]')}
  fils.poll.l <- setNames(lapply(organs,function(arg){paste(arg,'_',aux.l[[arg]],'.csv',sep='')}) , names(aux.l)) # filepaths [organs-list of filenames]
  for(c.organ in names(which(is.na(aux.l)))){fils.poll.l[[c.organ]] <- c()} # set to NA elements corresponding to organs without pure-pollutant profiles

 # Load pre-specified functions
  source('Modules/mollify.R')
  source('Modules/extrema.locs.R')
  source('Modules/dat.vs.pure.R')
  library(scales) # alpha transparency of plotting colors for scatter plots

# FIXED SETTINGS - DO NOT TOUCH UNLESS RESTRUCTURING CODEBASE -->>
####################################################################################################

   print('Settings complete')