#R code to calculate the life-cycle-NUE: 
library("Matrix")                       # Matrix calculations
library("MASS")                         # Inverted matrix

# BAsed on "Supply-and-Use framework" of Suuh and Yeh (2011)
# Based on source code received from Aimable Uwizeye (20161205)
# Modified and modularized for scenario calculations by Adrian Leip (20161208)


finputfactors<-function(n,s){
    #Total input relative to 1 final output unit
    i<-rev(sapply(1:n,function(x) Reduce("*",s[3:(4-x),1])))
    return(i)
}
fproducts<-function(inputfactors,product){
    #Total intput (scaled to original flows)
    totinputs<-c(inputfactors[2:nstages]*product,product)
    return(totinputs)
}
fflows<-function(nstages,stages,product){
    # Flows as in LCA
    flows<-fproducts(finputfactors(nstages,stages),product)*stages
    return(flows)
}
fnetout<-function(flows){
    # Net output = Total input per stage minus losses
    # .. Calculated from sum of useful outputs relative to products
    netout<-apply(flows[,3:5],1,sum)
    return(netout)
}
ftotinputs<-function(nstages,stages,product){
    #Total intput (scaled to original flows)
    totinputs<-finputfactors(nstages,stages)*product
    return(totinputs)
}
fstagesnue<-function(netout,totinput){
    nue<-netout/totinputs
    #nue = (82%, 70%, 72%)
    return(nue)
}
fnewinputs<-function(nstages,product,chain){
    mainflow<-matrix(rep(1,nstages*nstages),ncol=nstages,nrow=nstages)
    for(i in 1:(nstages-1)) mainflow[i,i+1]<-0
    #Input from recycling (scaled to original flows)
    rec<-product*chain
    recinputs<-rec*mainflow
    #New N inputs (from previous stage or added to stage)
    #.. exlcludes recycling from same of subsequent stages
    newinputs<-totinputs-apply(recinputs,2,sum)
    return(newinputs)
}
fextinputs<-function(nstages,chain,totinputs,product){
    mainflow<-matrix(rep(1,nstages*nstages),ncol=nstages,nrow=nstages)
    for(i in 1:(nstages-1)) mainflow[i,i+1]<-0
    #Input from recycling (scaled to original flows)
    rec<-product*chain
    recinputs<-rec*mainflow
    #New N inputs (from previous stage or added to stage)
    #.. exlcludes recycling from same of subsequent stages
    newinputs<-totinputs-apply(recinputs,2,sum)
    #External inputs - as new inputs but excluding input from previous stages
    extinputs<-totinputs-apply(rec,2,sum)
    return(extinputs)
}

lcanue<-function(nstages,product,chain,netout,extinputs){
    rec<-product*chain
    newresource<-extinputs%*%ginv(as.matrix(netout*Diagonal(nstages)-rec))
    print(extinputs)
    lcanue<-(1/newresource[1,nstages])*100
    #life_cycle_nue = 36%
    return(lcanue)
}


processflows<-c("input","loss","coproducts","residuals","stockchange","product")
nstages<-3
nflows<-6
lca<-matrix (rep(NA,nflows*nstages),nrow=nstages, ncol=nflows, byrow=F)
rec<-matrix (rep(0,nstages*nstages),nrow=nstages, ncol=nstages, byrow=F)

# Fill data from example supply chain
lca[1,]<-c(100+10+2+10,NA,40,10,0,50)
lca[2,]<-c(50,NA,5,10,0,20)
lca[3,]<-c(20,NA,0,2,0,12.5)

# Recycling flows (as shares of total co-products output of the stage)
rec[1,1]<-10 #/sum(lca[1,2:4],na.rm=TRUE)
rec[1,2]<-50 #/sum(lca[1,2:4],na.rm=TRUE)
rec[2,1]<-10 #/sum(lca[2,2:4],na.rm=TRUE)
rec[2,3]<-20 #/sum(lca[2,2:4],na.rm=TRUE)
rec[3,1]<-2  #/sum(lca[3,2:4],na.rm=TRUE)

# Losses are calculated as N not in any kind of products
lca[,2]<-lca[,1]-apply(lca[,c(3:nflows)],1,sum,na.rm=T)

# Stages indicates for each stage the flow relative to intended product
stages<-lca/lca[,nflows]
# Product is the quantity of nutrient in the intended product
# as of the initial supply chain description
product<-lca[nstages,nflows]
# Chain is the quantity of internal flows
# relative to the intended product (=1)
chain<-rec/product

chain1<-chain
stages1<-stages
source("resourceeffmatrixscen.r")

stagesori<-stages
chainori<-chain



#Increasing the share of recycling
stages1<-stagesori
chain1<-chainori
stages1[1,3]<-0.75
stages1[1,4]<-0.30
stages1[1,nflows]<-0.95
chain1[1,1]<-chain[1,1]*stages1[1,4]/stages[1,4]

source("resourceeffmatrixscen.r")

#Increasing the share of recycling
stages1<-stagesori
chain1<-chainori
stages1[1,3]<-0.7
stages1[1,4]<-0.40
stages1[1,nflows]<-0.9
chain1[1,1]<-chain[1,1]*stages1[1,4]/stages[1,4]
source("resourceeffmatrixscen.r")


#Increasing the share of co-products (assuming they have lower value)
#For example, part of the crop cannot be used as feed 
stages1<-stagesori
chain1<-chainori
stages1[1,3]<-1.0
stages1[1,nflows]<-0.8
source("resourceeffmatrixscen.r")

#Processing residuals are used as feed rather than as fertilizer
stages1<-stagesori
chain1<-chainori
chain1[3,1]<-0
chain1[3,2]<-0.16


n<-nstages
#First scale to output
chain1[,nstages]<-chain1[,nstages]/stages1[nstages,nflows]
stages1[nstages,]<-stages1[nstages,]/stages1[nstages,nflows]
mainflow<-matrix(rep(1,nstages*nstages),ncol=nstages,nrow=nstages)
for(i in 1:(nstages-1)) mainflow[i,i+1]<-0
#input from later stages
#this must be done recursively from nstages-1 until 2 (because 1st stage doesn't get input from previous ones)
rec1<-apply(chain1*mainflow,2,sum)
inputfactors1[nstages-1]<-(stages[nstages-1,1]*stages[nstages,1]-rec1[nstages-1])/stages[nstages,1]

for(j in 1:(nstages-1)) chain1[j,]<-chain1[j,]*inputfactors1[nstages-1]/stages[nstages-1,1]



