## Load all prelims and set organ
 setwd('/Users/herzog/Documents/Machine/WUR/CURR/IMMUNE/Part4/CODE')
 source('Modules/extrema.locs.R')
 source('Modules/settings.R')
 organ <- organs[1]
 source('Modules/presets.R')

## Compute start/end pure profiles
 source('Modules/pure.prof.tab.R') # outcome: pure.prof.l [(start,end)-list of (idx,x,H2O,D2O)-col df's]

## Identify global maxima of start/end pure profiles [(max,min)-list of (loc,idx) df's])
   print('pure pofiles')
 extr.pure.l <- list( # find extrema of start/end profiles
 	start=extrema.locs( data.frame(x=pure.prof.l$start$x,y=pure.prof.l$start$D2O) , smoothen=FALSE ) , # start
 	end=extrema.locs( data.frame(x=pure.prof.l$end$x,y=pure.prof.l$end$D2O) , smoothen=FALSE ) ) # end
 for(c.pos in c('start','end'))
 {
  aux.idx <- which.max(pure.prof.l[[c.pos]][extr.pure.l[[c.pos]]$max$idx,'D2O']) # global maximum index (in sub-list) [nat scalar]
  extr.pure.l[[c.pos]]$max <- list(loc = extr.pure.l[[c.pos]]$max$loc[aux.idx] , idx = extr.pure.l[[c.pos]]$max$idx[aux.idx] )
 }
   print('pure pofile final result')
   print(extr.pure.l)

## Plot pure profiles
 plot(t.g , pure.prof.l$start$D2O , typ='l' , col='blue' , # start profile
 	ylim=c(0,max(rbind(pure.prof.l$start$D2O,pure.prof.l$end$D2O))) )
 lines(extr.pure.l$start$max$loc , pure.prof.l$start[extr.pure.l$start$max$idx,'D2O'] , typ='p' , col='blue' , pch=19)
 lines(extr.pure.l$start$min$loc , pure.prof.l$start[extr.pure.l$start$min$idx,'D2O'] , typ='p' , col='blue' , pch=19)
 c.shft <- extr.pure.l$start$max$loc - extr.pure.l$end$max$loc # current shift (from maximum of reference pure profile) [multiple of dt.g]
 lines(t.g+c.shft , pure.prof.l$end$D2O , typ='l' , col='red') # end profile
 lines(extr.pure.l$end$max$loc+c.shft , pure.prof.l$end[extr.pure.l$end$max$idx,'D2O'] , typ='p' , col='red' , pch=19)
 lines(extr.pure.l$end$min$loc+c.shft , pure.prof.l$end[extr.pure.l$end$min$idx,'D2O'] , typ='p' , col='red' , pch=19)

for(c.n in setdiff(fils.v , c(ctrl.fils.df$fil,fils.poll.l[[organ]])))
{
   print(sprintf('dataset %s',c.n))
 # for(c.dupl in c(2,4))
 # {
   # print(sprintf('duplicate %f',c.dupl))
  aux.df <- dat.tab.l[[c.n]][,cols$D2O]
  aux.df <- sweep(aux.df , 2 , dt.g*colSums(aux.df) , '/')
  c.prof <- data.frame(x=t.g , y=rowMeans(aux.df)) # assign current data to c.prof [(x,y)-col df]
  c.prof$y <- c.prof$y/(dt.g*sum(c.prof$y)) # normalize profile to unit AUC [(x,y)-col df]
  c.extr.l <- extrema.locs(c.prof,smoothen=TRUE) # profile extrema [(max,min)-list of (loc,idx) lists])
  aux.idx <- which.min(abs(c.extr.l$max$loc-extr.pure.l$start$max$loc)) # extremum closest to reference point
  c.shft <- extr.pure.l$start$max$loc - c.extr.l$max$loc[aux.idx] # current shift (from maximum of reference pure profile) [multiple of dt.g]
    aux.clr <- 'gray'
    # aux.clr <- ifelse(c.dupl==2,'red','blue')
  lines(t.g+c.shft , mollify(c.prof)$y , typ='l' , col=aux.clr)
  lines(c.extr.l$max$loc+c.shft , mollify(c.prof)[c.extr.l$max$idx , 'y'] , typ='p' , col=aux.clr , pch=19)
  lines(c.extr.l$min$loc+c.shft , mollify(c.prof)[c.extr.l$min$idx , 'y'] , typ='p' , col=aux.clr , pch=19)
  readline(sprintf('%.02f',c.shft))
 # }

}