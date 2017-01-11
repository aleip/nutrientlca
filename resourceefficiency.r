doheader<-function(){
    library("Matrix")                       # Matrix calculations
    library("MASS")                         # Inverted matrix
}

supplychainsimple<-function(){
    processes<-c("Feed production",
                 "Livestock production",
                 "Food production")
    products<-c("feed","anim","food","feedco","animco","cropreso","cropresi","manure","foodres")
    resources<-c("N fertilizer")
    losses<-c("N losses")
    nenv<-length(resources)+length(losses)
    nproc<-length(processes)
    nprod<-length(products)
    
    #flagging co-products that are recycled
    recycling<-rep(0,length(products))
    recycling[6:9]<-1
    
    return(list(processes,nproc,products,nprod,resources,losses,nenv,recycling))
}

exampleaimable<-function(nproc,nprod,nenv){
    
    # Technology matrix
    A<-matrix(0,ncol=nproc,nrow=nprod)
    A[,1]<-c( 50,  0,   0,40,0,10,-10,-10,-2)
    # The vector is constructed such that consumption of recyling flows from other processes are not explicit to leave cross-dependencies flexible.
    # They are indicated with NA.
    # Required resources are coming from the resources (here mineral fertilizer)
    # This does NOT affect the calculated NUE per process or the lifecycle NUE based on material flow analysis.
    A[,1]<-c( 50,  0,   0,40,0,10,-10, NA,NA)
    A[,1]<-c( 25,  0,   0,20,0, 5, -5, NA,NA)
    A[,2]<-c(-50, 20,   0, 0,5, 0,  0, 10, 0)
    A[,3]<-c(  0,-20,12.5, 0,0, 0,  0,  0, 2)
    
    # Intervention matrix
    B<-matrix(0,ncol=nproc,nrow=nenv)
    B[,1]<-c(-100,22)
    B[,1]<-c(-112,22)
    B[,1]<-c(-56,11)
    B[,2]<-c(   0,15)
    B[,3]<-c(   0,5.5)
    
    # Demand
    f<-c(0,0,12.5) #Demand Vector required for main food chain
    
    return(list(A,B,f))
    
}

processmatrix<-function(techinter,supplychain){
    
    # Recover system definition
    processes<-supplychain[[1]]
    nproc<-supplychain[[2]]
    products<-supplychain[[3]]
    nprod<-supplychain[[4]]
    resources<-supplychain[[5]]
    losses<-supplychain[[6]]
    nenv<-supplychain[[7]]
    recycling<-supplychain[[8]]
    
    # Combine Technology and Intervention matrices to Process matrix 
    A<-as.data.frame(techinter[[1]])
    rownames(A)<-products
    colnames(A)<-processes
    B<-as.data.frame(techinter[[2]])
    rownames(B)<-interventions
    colnames(B)<-processes
    P<-rbind(A,B)
    colnames(P)<-processes
    
    return(list(A,B,P))
}

scaling<-function(A,f){
    # For the scaling factor only the main matrix counts
    # Co-products or recycling doesn't change the relative quantitive of the processes required
    # So far only P is scaled !
    s<-as.vector(ginv(as.matrix(A[1:ncol(A),]))%*%f)
    P[,]<-sapply(1:ncol(P),function(j) sapply(1:nrow(P),function(i) P[i,j]*s[j]))
    
    return(list(s,P))
}

fillrecycling<-function(recycling,P,nproc,nprod){
    temp<-apply(P[1:nprod,]*recycling,1,sum,na.rm=T)
    temp<-sapply(1:nproc,function(j) 
        sapply(1:nprod,function(i) 
            if(is.na(P[i,j])){-temp[i]}else{P[i,j]}))
    tempdiff<-as.vector(apply(temp,2,sum))-as.vector(apply(P[1:nrow(A),],2,sum,na.rm=T))
    
    P[1:nprod,]<-temp
    P[nprod+length(resources),]<-P[nprod+length(resources),]-tempdiff
    
    return(P)
}

flowanalysis<-function(P,nrod,recycling,goods,resources){
    chainflows<-sapply(1:nprod,function(x) if(round(sum(P[x,]),5)==0){which(P[x,]>0)>which(P[x,]<0)}else{FALSE})
    #flagging flows from the process they are generated (0 for imports)
    origin<-sapply(1:nprod,function(x) if(length(which(P[x,]>0))){which(P[x,]>0)}else{which(P[x,]<0)})
    target<-sapply(1:nprod,function(x) if(length(which(P[x,]<0))){which(P[x,]<0)}else{0})
    target[which(target==0&recycling==1)]<-origin[which(target==0&recycling==1)]
    exportflows<-target==0
    
    goodsproducts<-P[rownames(P)%in%goods,]
    goodsproducts[goodsproducts<0]<-0
    sumgoods<-apply(goodsproducts[rownames(P)%in%goods,],2,sum)
    goodstemp<-P[rownames(P)%in%goods,]
    goodstemp[goodstemp>0]<-0
    suminputs<--apply(goodstemp[rownames(P)%in%goods,],2,sum)
    suminputs<-suminputs-P[rownames(P)%in%resources,]
    sumresources<-P[rownames(P)%in%resources,]
    sumlosses<-P[rownames(P)%in%losses,]
    nue<-sumgoods/suminputs
    
    
    return(list(chainflows,origin,target,exportflows,sumlosses,sumresources,nue))
}


nutflowanalysi<-function(P,nproc,nprod,origin,resources){
    V<-matrix(rep(0,nproc**2),ncol=nproc,nrow=nproc)
    V<-t(sapply(1:nproc,function(x) apply(P[1:nprod,][origin==x,],2,sum)))
    r<-as.vector(as.matrix(-P[rownames(P)%in%resources,]))
    rintensity<-t(r)%*%ginv(V)
    lcanue<-1/rintensity[1,nproc]
    
    return(list(V,r,rintensity,lcanue))
}

allocationbyflow<-function(P,goods,origin){
    #Re-scale factors so that their sum is 1
    # -- note that for multiplication of columns by vector the number of rows must be the length of the vector
    goodsproducts<-P[rownames(P)%in%goods,]
    goodsproducts<-sapply(1:ncol(goodsproducts),function(x) goodsproducts[,x]*(origin==x))
    
    lambda<-t(t(goodsproducts)/apply(goodsproducts,2,sum))
    lambda<-t(t(goodsproducts))
    
    lambda<-t(t(lambda)/colSums(lambda))
    
    return(lambda)
}

allocationbyvalue1<-function(nproc,nprod){
    #Define relative 'value' (or whatever the basis for allocation is) over the goods produced in each process
    lam2<-matrix(0,ncol=nproc,nrow=nprod)
    # Main products 'double value' per unit of N
    lam2[1,1]<-P[1,1]*2
    lam2[2,2]<-P[2,2]*2
    lam2[3,3]<-P[3,3]*2
    # Exported products 'normal value
    lam2[4,1]<-P[4,1]*1
    lam2[5,2]<-P[5,2]*1
    # Recycled products backwards half value
    lam2[8,2]<-P[8,2]*0.5
    lam2[9,3]<-P[9,3]*0.5
    # Recycling to same process no burden
    lam2[6,1]<-P[6,1]*0
    
    lam2<-t(t(lam2)/colSums(lam2))
    
    # lam3 is the allocation to the products that are not flowing back into the supply chain
    # (thus exports or connecting flows)
    lam3<-lam2*(!chainflows)
    lam3<-t(t(lam3)/colSums(lam3))
    
    return(list(lam2,lam3))
    
}

calcburden<-function(P,nproc,nprod,lambda,lambdamain,chainflows,target){
    
    # Burden must be total losses!
    # Distribute the direct burden (total losses) of each process 
    # over the products of the process including those that are recycled
    rburd<-which(interventions%in%losses)+nprod
    rreso<-which(interventions%in%resources)+nprod
    dirburden<-sapply(1:nproc,function(j) 
        sapply(1:(nprod),function(i)-lambda[i,j]*as.vector(P[rburd,j])))
    
    # Resources must be calculated as N in product + N in losses
    # Calculate the resources each product requires to match
    # the equation: Resources = Products + Losses
    # Note - in the example there is no stock changes so far,
    #        equations would need to be adapted as stock changes 
    #        are resources but count as goods without burden
    A<-P[1:nprod,]
    dirresour<-dirburden
    lcells<-dirresour<0
    dirresour[lcells]<-A[lcells]-dirburden[lcells]
    
    # Calculate the embedded burden in recycling flows
    # This formula calculates how much burden of a chainflow needs 
    # to be added to the target process
    recburden<-colSums(sapply(1:nproc,function(x) rowSums(dirburden*chainflows*(target==x))))
    
    # The embedded burden is distributed ...
    embburden<-sapply(1:nproc,function(j) 
        sapply(1:(nprod),function(i) lambda[i,j]*recburden[j]))
    
    # ... and added to the direct burden
    # This yields the burden for each process
    # Through recycling the burden of processes that receive recycling flows
    # is increased and those recycling is decreased
    burdenprocesses<-dirburden+embburden
    burdenprocesses[chainflows,]<-0
    
    # Finally the burden is carried along the chain through
    # the 'connecting' flows
    # The embedded burden from previous processes is distributed
    # amongst the exported goods and the flows to successive processes
    
    # Attention! this is programmed only for a three-process system
    # Need generalization !!
    burdenproducts<-burdenprocesses
    burdenproducts[,1]<-burdenprocesses[,1]
    burdenproducts[,2]<-burdenprocesses[,2]+burdenproducts[1,1]*lambdamain[,2]
    burdenproducts[,3]<-burdenprocesses[,3]+burdenproducts[2,2]*lambdamain[,3]
    burdenproducts[1,1]<-0
    burdenproducts[2,2]<-0
    
    # Resources must be calculated as N in product + N in losses
    resourcesproducts<-burdenproducts
    lcells<-round(burdenproducts,5)<0
    resourcesproducts[lcells]<-A[lcells]-burdenproducts[lcells]
    finproducts<-A
    finproducts[!lcells]<-0
    
    
    # Nutrient Use Efficiency
    nueproducts<-finproducts/resourcesproducts
    
    # Footprints
    lossfactors<--burdenproducts/finproducts
    
    return(list(dirburden,recburden,embburden,
                burdenproducts,resourcesproducts,finproducts,nueproducts,
                lossfactors))
}


#Matrix of flows which are 'too many' (not part of the main food chain)
# ... from Heijungs adn Suh, but here probably not needed
# Calculation requires a square Technology matrix
# Fill missing columns on the basis of existing data
Afull<-function(A){
    n<-ncol(A)
    m<-nrow(A)
    grows2add<-m-n
    g<-matrix(rep(0,ncol=(m-n)*(m-n)),ncol=(m-n),nrow=(m-n))
    gtemp<-matrix(rep(0,ncol=(m-n)*n),ncol=(m-n),nrow=n)
    fsum<-apply(A,1,sum)[(n+1):m]
    for(i in 1:length(fsum)){g[i,i]<--fsum[i]}
    g<-rbind(gtemp,g)
    g<-cbind(A,g)
    return(g)
}
fplus<-function(A,f){
    g<-c(f,rep(0,nrow(A)-ncol(A)))
    return(g)
}
