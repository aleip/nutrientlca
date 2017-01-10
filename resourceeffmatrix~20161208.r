#R code to calculate the life-cycle-NUE: 
library("Matrix")                       # Matrix calculations
library("MASS")                         # Inverted matrix

lcnue<-function(lca,rec){
    

    # Net output = Total input per stage minus losses
    netout<-lca[,1]+lca[,4]-lca[,5]
    
    nue<-100*netout/lca[,1]
    #nue = (82%, 70%, 72%)
    
    newinput<-lca[,1]-apply(rec,2,sum)
    newresource<-newinput%*%ginv(as.matrix(netout*Diagonal(nstages)-rec))
    lcanue<-(1/newresource[1,3])*100
    #life_cycle_nue = 36%
}

stagesnue<-function(products,stages,chain){
    
    
    #products[2]<-stages[3,1]*products[3]/stages[3,2]
    #products[1]<-stages[2,1]*products[2]/stages[2,1]
    
    # Net output = Total input per stage minus losses
    netout<-products/stages[,2]
    #[1] 100.0  35.0  14.5
    
    mainflow<-matrix(rep(1,nstages*nstages),ncol=nstages,nrow=nstages)
    for(i in 1:(nstages-1)) mainflow[i,i+1]<-0
    #Total intput (scaled to original flows)
    totinputs<-netout*stages[,1]
    #[1] 122  50  20
    
    #Input from recycling (scaled to original flows)
    recinputs<-products[3]*chain*mainflow
    #New N inputs (from previous stage or added to stage)
    #.. exlcludes recycling from same of subsequent stages
    newinputs<-netout*stages[,1]-apply(recinputs,2,sum)
    #External inputs - as new inputs but excluding input from previous stages
    extinputs<-lca[,1]-apply(rec,2,sum)
    
    nue<-netout/totinputs
    #nue = (82%, 70%, 72%)
    
    newresource<-extinputs%*%ginv(as.matrix(netout*Diagonal(nstages)-rec))
    lcanue<-(1/newresource[1,3])*100
    #life_cycle_nue = 36%
}


processflows<-c("input","product","coproducts","stockchange","loss")
nstages<-3
nflows<-5
lca<-matrix (rep(NA,nflows*nstages),nrow=nstages, ncol=nflows, byrow=F)
rec<-matrix (rep(0,nstages*nstages),nrow=nstages, ncol=nstages, byrow=F)

# Fill data from example supply chain
lca[1,]<-c(100+10+2+10,50,40+10,0,NA)
lca[2,]<-c(50,20,5+10,0,NA)
lca[3,]<-c(20,12.5,2,0,NA)

# Recycling flows (as shares of total co-products output of the stage)
rec[1,1]<-10 #/sum(lca[1,2:4],na.rm=TRUE)
rec[1,2]<-50 #/sum(lca[1,2:4],na.rm=TRUE)
rec[2,1]<-10 #/sum(lca[2,2:4],na.rm=TRUE)
rec[2,3]<-20 #/sum(lca[2,2:4],na.rm=TRUE)
rec[3,1]<-2  #/sum(lca[3,2:4],na.rm=TRUE)

# Losses are calculated as N not in any kind of products
lca[,5]<-lca[,1]-apply(lca[,2:4],1,sum,na.rm=T)

# Stages indicates for each stage the flow relative to the sum
# of useful products (product, co-products + residuals + recycling, stock changes)
stages<-lca/apply(lca[,2:4],1,sum,na.rm=T)
# Product is the quantity of nutrient in the intended product
# as of the initial supply chain description
products<-lca[,2]
# Chain is the quantity of internal flows
# relative to the intended product (=1)
chain<-rec/products[3]



mylca<-lcnue(lca,rec)
mylca<-stagesnue(products = products,stages = stages,chain = chain)

# Increasing the share of co-products (assuming they have lower value)
#stages[1,3]<-0.75
#stages[1,2]<-1-stages[1,3]

