source("resourceefficiency_header.r")
# A. Process example 
S<-supplychainsimple()
E<-exampleaimable(S$nproc,S$nprod,S$nenv)
#E<-processmatrix(E,S)
E$P<-scaling(E$A,E$B,E$f)
E$P<-fillrecycling(S$recycled,S$resources,E$P,S$nproc,S$nprod)
flows<-flowanalysis(E$P,S$nprod,S$recycled,S$goods,S$resources,S$losses)

# B. Calculation of material flow analysis
nutflow<-nutflowanalysi(E$P,S$nproc,S$nprod,flows$origin,S$resources,S$losses,S$prows)

# C. Calculation of efficiency acc to allocation
lambda<-allocationbyflow(E$P,S$goods,flows$origin,S$recycled,S$prows) 
flowmatrix<-f_flowmatrix(S$nproc,S$nprod,lambda,flows$origin,flows$target)
flowfin<-f_flowfin(flowmatrix)
flowrec<-f_flowrec(flowmatrix)
recfate<-f_recfate(S$nproc,flowfin,flowrec)

burden<-calcburden(E$P,S$nproc,S$nprod,lambda,recfate,flows$chainflows,flows$target,flows$origin,flows$exportflows,S$resources,S$losses)

print(1/nutflow$rintensity)
print(colSums(burden$finproducts)*nutflow$rintensity)
#print(nutflow$lcanue)
print(colSums(burden$nueproducts))
print(colSums(burden$finproducts)/colSums(burden$nueproducts))


