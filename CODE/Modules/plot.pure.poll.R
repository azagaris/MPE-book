####################################################################################################
# 'plot.pure.poll.profs.R' [Antonios Zagaris, Wageningen Bioveterinary Research (WUR) 2017]
#
# This routine plots the pure and pure-pollutant profiles for the organ currently in the workspace
# (and set by main.R). The profiles are normalized to have unit AUC over the rough interval (cf.
# settings.R) and are only plotted over it. Also plotted are the cutoffs for the start/end pure
# profiles, i.e. the lower and upper bounds (as well as the valley) for those profiles and for the
# current organ. Although these cutoffs are included (as vertical lines), all portrayed profiles are
# normalized over the rough interval and not over the corresponding interval bounded by the cutoffs.
#
#   INPUTS
# None.
#
#   OUTPUTS
# * Plots of all (H2O/D2O) datasets (blue/red) for the current organ, together with (H2O/D2O) pure
#   profiles (light blue/pink), the pollutant profile (green) and vertical lines at the (sliding;  
#   per dataset) left/right (H2O/D2O) cutoffs and valley (blue/red)
#
# * Various plot settings, used also futher down the pipeline
#
#   DEPENDENCIES
# <- settings
# <- presets
# <- all modules before it in main.R
# -> all modules after it in main.R
####################################################################################################
# PLOT SETTINGS -->>

 plot.eps.df <- data.frame(x=2E-3,y=2E-2) # plotting buffers [single row, (x,y)-col df]
 plot.max.y <- (1+plot.eps.df$y)*max(pure.prof.l$start[,c('H2O','D2O')] , pure.prof.l$end[,c('H2O','D2O')]) # max y [scalar]
 plot.label.lims <- list( # start/end coordinates for inserting 'start', 'end' plot labels [(start,end)-list of single-row (x,y)-df's]
   start = data.frame( x=t.g[which.max(pure.prof.l$start$H2O)] , y=(1+plot.eps.df$y)*max(pure.prof.l$start$H2O) ) ,
   end = data.frame( x=t.g[which.max(pure.prof.l$end$H2O)] , y=(1+plot.eps.df$y)*max(pure.prof.l$end$H2O) ) )
 plot.label.lims$start$x <- (1 + 2*((plot.label.lims$start$x>plot.label.lims$end$x) - 1/2)*plot.eps.df$x)*plot.label.lims$start$x # 'start' label position
 plot.label.lims$end$x <- (1 - 2*((plot.label.lims$start$x>plot.label.lims$end$x) - 1/2)*plot.eps.df$x)*plot.label.lims$end$x # 'end' label position

# PLOT SETTINGS <<--
####################################################################################################
# PLOTS -->>

if(plot.swtc){ # check plot switch (only (re)plot if switch is on)
   # pdf(paste('../Results/',organ,'/pure_and_poll_profiles_',organ,'.pdf',sep='')) # open plotting connection

 plot(t.g , pure.prof.l$start$H2O , typ='l' , col='blue' ,
      ylim=c(0,plot.max.y) , xlab='time' , ylab='density' , main=organ) # pure start H2O profile
 text(plot.label.lims$start$x , plot.label.lims$start$y , 'start' , col='black')
 text(plot.label.lims$end$x , plot.label.lims$end$y , 'end' , col='grey')
 lines(t.g , pure.prof.l$end$H2O , typ='l' , col='lightblue') # pure end H2O profile
 lines(t.g , pure.prof.l$start$D2O , typ='l' , col='red') # pure start D2O profile
 lines(t.g , pure.prof.l$end$D2O , typ='l' , col='pink') # pure end D2O profile

 lines(t.g , poll.prof.df$y , typ='l' , col='green') # pollutant profile

   # dev.off() # close plotting connection
} # end of check plot switch

# PLOTS <<--
####################################################################################################