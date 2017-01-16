#E<-processmatrix(E,S)
s<-scaling(E$A,E$B,E$f)
E$P<-s$P
E$P<-fillrecycling(E,S)
flows<-flowanalysis(E$P,S$nprod,S$recycled,S$goods,S$resources,S$losses)
# B. Calculation of material flow analysis
nutflow<-nutflowanalysi(E$P,S$nproc,S$nprod,flows$origin,S$resources,S$losses,S$prows)

# C. Calculation of efficiency acc to allocation
lambda<-allocationbyflow(E$P,S$goods,flows$origin,S$recycled,S$prows) 
flowmatrix<-f_flowmatrix(S$nproc,S$nprod,lambda,flows$origin,flows$target)
flowfin<-f_flowfin(flowmatrix)
flowrec<-f_flowrec(flowmatrix)
chainfractions<-f_chainfraction(S$nproc,flowfin)
recfate<-f_recfate(S$nproc,flowrec,chainfractions)
burden<-calcburden(E,S,flows,lambda,recfate)
