#E<-processmatrix(E,S)
s<-scaling(E$A,E$B,E$f)
E$P<-s$P
E$P<-fillrecycling(E,S)
flows<-flowanalysis(E,S)
# B. Calculation of material flow analysis
nutflow<-nutflowanalysi(E,S,flows)

# C. Calculation of efficiency acc to allocation
if(dolambda=="byflow"){
    lambda<-allocationbyflow(E,S,flows) 
}else if(dolambda=="byvalue1"){
    lambda<-allocationbyvalue1(E$P,S$nprod,S$nproc) 
}else if(dolambda=="byvalue2"){
    lambda<-allocationbyvalue2(E$P,S$nprod,S$nproc) 
}
flowmatrix<-f_flowmatrix(S$nproc,S$nprod,lambda,flows$origin,flows$target)
flowfin<-f_flowfin(flowmatrix)
flowrec<-f_flowrec(flowmatrix)
chainfractions<-f_chainfraction(S$nproc,flowfin)
recfate<-f_recfate(S,flowrec,chainfractions)
burden<-calcburden(E,S,flows,lambda,recfate)
