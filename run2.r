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

lambda<-allocationbyflow(E$P,S$goods,flows$origin,S$recycled,S$prows) 
flowmatrix<-f_flowmatrix(S$nproc,S$nprod,lambda,flows$origin,flows$target)
flowfin<-f_flowfin(flowmatrix)
flowrec<-f_flowrec(flowmatrix)
recfate<-f_recfate(E$P,S$nproc,S$nprod,lambda,flows$chainflows,flows$exportflows)

burden<-calcburden(E$P,S$nproc,S$nprod,lambda,recfate,flows$chainflows,flows$target,flows$origin,flows$exportflows,S$resources,S$losses)

burden$nueproducts
