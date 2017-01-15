library("Matrix")                       # Matrix calculations
library("MASS")                         # Inverted matrix

rm(list=objects())

source("genfunctions.r")
source("supplychaindefs.r")
source("supplychainexamples.r")
source("resourceefficiency.r")

S<-supplychainsimple()
E<-exampleaimable(S$nproc,S$nprod,S$nenv)
#E<-processmatrix(E,S)
E$P<-scaling(E$A,E$B,E$f)
E$P<-fillrecycling(S$recycled,S$resources,E$P,S$nproc,S$nprod)
flows<-flowanalysis(E$P,S$nprod,S$recycled,S$goods,S$resources,S$losses)
nutflow<-nutflowanalysi(E$P,S$nproc,S$nprod,flows$origin,S$resources,S$losses,S$prows)

nutflow$lcanue

