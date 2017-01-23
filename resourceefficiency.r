processmatrix<-function(E,S){
    # Recover system definition
    # Combine Technology and Intervention matrices to Process matrix 
    A<-as.data.frame(E$A)
    rownames(A)<-S$products
    colnames(A)<-S$processes
    B<-as.data.frame(E$B)
    rownames(B)<-c(S$resources,S$losses)
    colnames(B)<-S$processes
    P<-rbind(A,B)
    colnames(P)<-S$processes
    
    return(list(A=A,B=A,P=P))
}
scaling<-function(A,B,f){
    # For the scaling factor only the main matrix counts
    # Co-products or recycling doesn't change the relative quantitive of the processes required
    # So far only P is scaled !
    P<-rbind(A,B)
    s<-as.vector(ginv(as.matrix(A[1:ncol(A),]))%*%f)
    P[,]<-sapply(1:ncol(P),function(j) sapply(1:nrow(P),function(i) P[i,j]*s[j]))
    
    return(list(s=s,P=P))
}

fillrecycling<-function(E,S){
    recycled<-S$recycled
    resources<-S$resources
    P<-E$P
    nproc<-S$nproc
    nprod<-S$nprod
    lrow<-recycled==1
    lrow<-lrow|is.na(apply(P[1:nprod,],1,sum))
    temp<-apply(P[1:nprod,]*lrow,1,sum,na.rm=T)
    temp<-sapply(1:nproc,function(j) 
        sapply(1:nprod,function(i) 
            if(is.na(P[i,j])){-temp[i]}else{P[i,j]}))
    tempdiff<-as.vector(apply(temp,2,sum))-as.vector(apply(P[1:nprod,],2,sum,na.rm=T))
    
    P[1:nprod,]<-temp
    P[nprod+length(resources),]<-P[nprod+length(resources),]-tempdiff
    return(P)
}

flowanalysis<-function(E,S){
    #E$P,S$nprod,S$recycled,S$goods,S$resources,S$losses
    #P,nprod,recycled,goods,resources,losses
    P<-E$P
    nprod<-S$nprod
    recycled<-S$recycled
    products<-S$products
    goods<-S$goods
    resources<-S$resources
    losses<-S$losses
    chainflows<-sapply(1:nprod,function(x) if(round(sum(P[x,]),5)==0 & length(which(P[x,]!=0)>0)){which(P[x,]>0)>which(P[x,]<0)}else{FALSE})
    
    #flagging flows from the process they are generated (0 for imports)
    origin<-sapply(1:nprod,function(x) if(length(which(P[x,]>0))){which(P[x,]>0)}else if(length(which(P[x,]<0))){which(P[x,]<0)}else{0})
    target<-sapply(1:nprod,function(x) if(length(which(P[x,]<0))){which(P[x,]<0)}else{0})
    target[which(target==0&recycled==1)]<-origin[which(target==0&recycled==1)]
    exportflows<-target==0
    
    goodsproducts<-P[rownames(P)%in%goods,]
    goodsproducts[goodsproducts<0]<-0
    sumgoods<-apply(goodsproducts,2,sum)
    goodstemp<-P[rownames(P)%in%goods,]
    goodstemp[goodstemp>0]<-0
    suminputs<--apply(goodstemp,2,sum)
    suminputs<-suminputs-P[rownames(P)%in%resources,]
    sumresources<-P[rownames(P)%in%resources,]
    sumlosses<-P[rownames(P)%in%losses,]
    nue<-sumgoods/suminputs
    
    names(chainflows)<-products
    names(origin)<-products
    names(target)<-products
    
    return(list(chainflows=chainflows,
                origin=origin,
                target=target,
                suminputs=suminputs,
                sumgoods=sumgoods,
                exportflows=exportflows,
                sumlosses=sumlosses,
                sumresources=sumresources,
                nue=nue))
}

nutflowanalysi<-function(E,S,flows){
    P<-E$P
    nproc<-S$nproc
    nprod<-S$nprod
    origin<-flows$origin
    resources<-S$resources
    losses<-S$losses
    waste<-S$waste
    goods<-S$goods
    prows<-S$prows
    
    sel<-which(rownames(P)%in%goods)
    #print(apply(as.matrix(P[sel,][origin[sel]==3,]),2,sum))
    #print(P[sel,][origin[sel]==3,,drop=FALSE])
    V<-matrix(rep(0,nproc**2),ncol=nproc,nrow=nproc)
    V<-t(sapply(1:nproc,function(x) apply(P[sel,][origin[sel]==x,,drop=FALSE],2,sum)))
    r<-as.vector(as.matrix(-P[prows%in%resources,]))
    e<-as.vector(as.matrix(-P[prows%in%losses,]))
    rintensity<-t(r)%*%ginv(V)
    eintensity<-t(e)%*%ginv(V)
    lcanue<-1/rintensity[1,nproc]
    
    rownames(V)<-colnames(P)
    colnames(rintensity)<-colnames(P)
    colnames(eintensity)<-colnames(P)
    
    return(list(V=V,r=r,rintensity=rintensity,lcanue=lcanue))
}



allocationbyflow<-function(E,S,flows){
    #Re-scale factors so that their sum is 1
    # -- note that for multiplication of columns by vector the number of rows must be the length of the vector
    P<-E$P
    goods<-S$goods
    origin<-flows$origin
    recycled<-S$recycled
    prows<-S$prows
    products<-S$products
    nprod<-S$nprod
    
    sel<-which(products%in%goods)
    sel<-1:nprod
    goodsproducts<-P[sel,]
    goodsproducts<-sapply(1:ncol(goodsproducts),function(x) goodsproducts[,x]*(origin[sel]==x))
    
    lam2<-t(t(goodsproducts)/apply(goodsproducts,2,sum))
    lam2<-t(t(goodsproducts))
    lam2[recycled[sel]==1,]<-0
    
    lam2<-t(t(lam2)/colSums(lam2))
    colnames(lam2)<-colnames(P)
    return(lambda=lam2)
}

f_flowmatrix<-function(nproc,nrpod,lambda,origin,target,goods){
    flowmatrix<-matrix(0,ncol=nproc,nrow=nproc)
    flowmatrix<-sapply(1:nproc,function(j) sapply(1:nproc,function(i)
        if(i!=j){
            sum(lambda[origin==i&target==j,])
        }else{
            sum(lambda[origin==i&target==0,])
        }
    ))
    #sel<-which(rownames(lambda)%in%goods)
    #print(apply(as.matrix(P[sel,][origin[sel]==3,]),2,sum))
    #print(P[sel,][origin[sel]==3,,drop=FALSE])
    #V<-matrix(rep(0,nproc**2),ncol=nproc,nrow=nproc)
    #V<-t(sapply(1:nproc,function(x) apply(lambda[sel,][origin[sel]==x,,drop=FALSE],2,sum)))
    
    return(flowmatrix)
}
f_flowfin<-function(flowmatrix){
    flowfin<-flowmatrix*upperTriangle(flowmatrix)
    #flowfin<-(flowfin/rowSums(flowfin))
    return(flowfin)
}
f_flowrec<-function(flowmatrix){
    #share recycled are not recycled    
    flowrec<-flowmatrix*(1-upperTriangle(flowmatrix))
    #for(i in 1:ncol(flowrec)) flowrec[i,i]=1-sum(flowrec[i,])
    return(flowrec)
}

f_chainfraction<-function(nproc,flowfin){
    #nproc<-S$nproc
    chainfractions<-matrix(0,ncol=nproc,nrow=nproc)
    
    # Approach: separate 'direct' chain-flow ignoring recycling flows and
    # then calculate recycling flows.
    
    # First calculate the distribution of burden that are transported
    # along the supply chain.
    # This loop calculates 'chainfractions' of burden that are 'caused' by
    # a certain process i (rows) over all processes j (columns)
    # The burden is attributed to the sum of products that leave each respective
    # process.
    
    for(j in nproc:1){
        # Loop over all last-but-one origins for the target
        for(i in j:1){
            # The fraction arriving at target is the direct flow ... 
            chainfractions[i,j]<-flowfin[i,j]*flowfin[j,j]
            if(j==i)chainfractions[i,j]<-flowfin[j,j]
            
            if(i<j-1){
                for(ii in (i+1):(j-1)){
                    #chainfractions[ii,i]<-chainfractions[ii,i]+flowfin[ii,i]*if(i<j){chainfractions[i,j]}else{1}
                    
                    chainfractions[i,j]<-chainfractions[i,j] +
                        # ... plus the fraction from origin-to-midorigin*fraction from 
                        # mid-origin(already calculated) to target
                        flowfin[i,ii]*chainfractions[ii,j]
                }
            }
        }
    }
    chainfractions<-chainfractions/rowSums(chainfractions)
    
    
    
    
    return(chainfractions)
}
f_recfate<-function(S,flowrec,chainfractions){
    nproc<-S$nproc
    #nproc<-S$nproc
    finfractions<-matrix(0,ncol=nproc,nrow=nproc)
    for(i in 1:nproc){
        
        for(ii in 1:max(1,(i-1))){
            finfractions[i,]<-finfractions[i,]+
                
                # The recycled part enters in a previous step
                flowrec[i,ii]*finfractions[ii,]
        }
        
        finfractions[i,]<-finfractions[i,]+
            # The non-recycled part is distributed as in defined in  
            # chainfractions along the chain
            (1-sum(flowrec[i,]))*chainfractions[i,]
    }
    #print(flowrec)
    #print(chainfractions)
    #print(finfractions)
    #return(recfate=finfractions)
    return(recfate=chainfractions)
}

calcmatrix<-function(V,lambda,dirburden){
    #Retrieve the sum of net goods per process
    sumnetgoods<-colSums(V*diag(3))
    temp<-t(t(lambda)*sumnetgoods)
    
    #Construct matrix
    E<-V
    E[,]<-0
    E[1,2]<--dirburden[1,1]
    E[2,1]<--dirburden[8,2]
    E[3,1]<--dirburden[9,3]
    E[2,3]<--dirburden[2,2]
    E<-E-diag(3)*colSums(dirburden)
    
    lossintensity<-t(sumnetgoods)%*%ginv(E)
    
    lossintensity%*%E
}


calcburden<-function(E,S,flows,lambda,recfate){
    P<-E$P
    nproc<-S$nproc
    nprod<-S$nprod
    chainflows<-flows$chainflows
    target<-flows$target
    origin<-flows$origin
    exportflows<-flows$exportflows
    resources<-S$resources
    losses<-S$losses
    # Burden must be total losses!
    # Distribute the direct burden (total losses) of each process 
    # over the products of the process including those that are recycled
    interventions<-c(resources,losses)
    rburd<-which(interventions%in%losses)+nprod
    rreso<-which(interventions%in%resources)+nprod
    temp<-P[1:nprod,]
    temp[,]<-0
    temp2<-colSums(temp)
    dirburden<-temp
    dirburden[,]<-sapply(1:nproc,function(j) 
        sapply(1:(nprod),function(i)lambda[i,j]*as.vector(P[rburd,j])))
    
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
    # Add burden from original chainflows
    #recburden<-recburden+c(0,colSums(dirburden[1:(nproc-1),])[1:(nproc-1)])
    recburden<-recburden+colSums(dirburden*(!chainflows))
    # The embedded burden is distributed over the processes...
    embburden<-recburden*recfate
    embburdenprc<-colSums(embburden)
    # Distribute over exported products
    exproducts<-A*exportflows
    #print(t(embburdenprc*t(exproducts)))
    embburdenprd<-t(embburdenprc*t(exproducts))
    embburdenprd<-t(t(embburdenprd)/colSums(exproducts))
    
    # ... and added to the direct burden
    # This yields the burden for each process
    # Through recycling the burden of processes that receive recycling flows
    # is increased and those recycling is decreased
    #burdenproducts<-dirburden*exportflows+embburdenprd
    burdenproducts<-embburdenprd
    # Resources must be calculated as N in product + N in losses
    finproducts<-A*exportflows
    resourcesproducts<-finproducts+burdenproducts
    nfinprod<-sum(finproducts>0)
    cfinprod<-which(finproducts>0,arr.ind=TRUE)
    
    # Nutrient Use Efficiency
    nueproducts<-finproducts*0
    lcell<-finproducts!=0
    nueproducts[lcell]<-finproducts[lcell]/resourcesproducts[lcell]
    
    # Footprints
    lossfactors<-finproducts*0
    lossfactors[lcell]<--burdenproducts[lcell]/finproducts[lcell]
    
    return(list(dirburden=dirburden,
                dirresour=dirresour,
                recburden=recburden,
                embburden=embburden,
                embburdenprd=embburdenprd,
                burdenproducts=burdenproducts,
                resourcesproducts=resourcesproducts,
                finproducts=finproducts,
                nueproducts=nueproducts,
                lossfactors=lossfactors,
                nfinprod=nfinprod,
                cfinprod=cfinprod))
}

f_reffanalysis<-function(
    supplydef="default",
    supplyexe="aimable",
    supplyall="NA"
){
    # Function to generate a Lists (examp) of Lists (burden, nutflow, ...) 
    # with the results.
    # The number of scenarios run is given by the length of the vector of 
    # 'examples' with flow rates (supplyexe).
    # The other two vectors supplydef and supplyall are indicating the 
    # system definitions and allocation rules to be used. Those vectors can
    # be shorter (or just contain one value) in which case the last value is
    # used for all further scenarios
    # 
    # Allocations:
    # NA for applying MFA aprroach
    # 'byflow','byvalue1','byvalue2' etc for LCA appraoch
    
    nexamples<-length(supplyexe)
    
    for(i in 1:nexamples){
        
        stype<-if(length(supplydef)<i){supplydef[length(supplydef)]}else{supplydef[i]}
        etype<-supplyexe[i]
        dolambda<-if(length(supplyall)<i){supplyall[length(supplyall)]}else{
            supplyall[i]}
        if(!is.na(dolambda)) {if(dolambda=="NA") {dolambda<-NA}}
        
        cat(i,stype,etype,dolambda,"\n")
        S<-supplychainsimple(stype)
        E<-supplyvalues(etype,S)
        
        s<-scaling(E$A,E$B,E$f)
        E$P<-s$P
        E$P<-fillrecycling(E,S)
        flows<-flowanalysis(E,S)
        # B. Calculation of material flow analysis
        if(is.na(dolambda)){
            nutflow<-nutflowanalysi(E,S,flows)
            curex<-list(stype=stype,etype=etype,dolambda=dolambda,
                        S=S,E=E,flows=flows,nutflow=nutflow)
        }else{
            # C. Calculation of efficiency acc to allocation
            if(dolambda=="byflow"){
                lambda<-allocationbyflow(E,S,flows) 
            }else if(dolambda=="byvalue1"){
                lambda<-allocationbyvalue1(E$P,S$nprod,S$nproc) 
            }else if(dolambda=="byvalue2"){
                lambda<-allocationbyvalue2(E$P,S$nprod,S$nproc) 
            }
            #print(lambda)
            flowmatrix<-f_flowmatrix(S$nproc,S$nprod,lambda,flows$origin,flows$target)
            #print(flowmatrix)
            flowfin<-f_flowfin(flowmatrix)
            print("flowfin")
            print(flowfin)
            flowrec<-f_flowrec(flowmatrix)
            #print('flowrec')
            #print(flowrec)
            chainfractions<-f_chainfraction(S$nproc,flowfin)
            
            print("chainfractions")
            print(chainfractions)
            recfate<-f_recfate(S,flowrec,chainfractions)
            print(recfate)
            burden<-calcburden(E,S,flows,lambda,recfate)
            curex<-list(stype=stype,etype=etype,dolambda=dolambda,
                        S=S,E=E,flows=flows,burden=burden,
                        lambda=lambda,flowmatrix=flowmatrix,flowfin=flowfin,
                        flowrec=flowrec,chainfractions=chainfractions,recfate=recfate)
        }
        if(i==1){examp<-curex}else if(i==2){examp<-list(examp,curex)}else{
            examp[[i]]<-curex
        }
    }
    return(examp)
}

comparison<-function(a,what2,which2=NULL){
    if(is.null(which2))which2<-1:length(a)
    if(what2%in%"P"){
        for(i in which2){
            b<-a[[i]]
            cat("\nProcess matrix",i,"\n")
            print(b$E$P)
        }
    }
    if(what2%in%"nue"){
        cat("\nNutrient Use Efficiency - system definition - example\n")
        for(i in which2){
            b<-a[[i]]
            if(is.na(b$dolambda)){
                approach<-"MFA"
                cat(i,approach,round(as.vector(1/b$nutflow$rintensity),2),
                    "-",b$stype,"-",b$etype,"\n")
            }else{
                approach<-"LCA"
                b$burden$nueproducts[b$burden$nueproducts==0]<-NA
                cat(i,approach,round(apply(b$burden$nueproducts,2,mean,na.rm=TRUE),5),
                    "-",b$stype,"-",b$etype,"-",b$dolambda,"\n")
            }
        }
    }
}

describe<-function(X){
    n<-which(X!=0,arr.ind = TRUE)
    m<-data.frame(matrix(0,nrow=nrow(n),ncol=3))
    colnames(m)<-c("Process","Product","Value")
    m[,2]<-rownames(n)
    m[,1]<-colnames(X)[n[,2]]
    m[,3]<-sapply(1:nrow(m),function(x) X[n[x,1],n[x,2]])
    return(m)
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
temploop<-function(V){
    
    i=1
    emb<-0
    for(x in 1:1000){
        
        a=50/90
        a2=40/90
        b=10/35
        b2=5/35
        c=20/35
        d=2/14.5
        e=12.5/14.5
        
        
        # The recycled input is returned through recycled outputs
        retur=(a*b+a*c*d)
        # All must be leaving the system in the final products
        leave=(a2+a*b2+a*c*e)
        retur2<-retur*retur
        leave2<-leave+retur*leave
        #
        retur3<-retur**3
        leave3<-leave+retur*leave+retur**2*leave
        
        F1<-a2/leave
        F2<-a*b2/leave
        F3<-a*c*e/leave
        #Ft must give 1
        Ft<-F1+F2+F3
        
        exitprod<-lambda*(target==0)
        cmult<-diag(nproc)
        cmult[2:nproc,2:nproc]<-lambda[1:(nproc-1),1:(nproc-1)]
        cmult<-colSums(cmult)
        
        
        j=a*b+a*c*d
        o=a*b+a*c*e
        emb<-emb+j
        if(j<0.00001) break
    }
}
