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
f_recfatelambda<-function(E,S,flows,lambda,flowfin){
    P<-E$P
    chainflows<-flows$chainflows
    nproc<-S$nproc
    nprod<-S$nprod
    
    lam3<-lambda*(!(chainflows))
    lam3<-t(t(lam3)/colSums(lam3))
    
    
    #Calculate the destination of an input to process 1
    #that is not recycled
    
    # Create a matrix that 'chain fractions', ie. the 
    # share of process j (column) 
    # that is arriving at process i (rows)
    chainfractions<-matrix(0,ncol=nproc,nrow=nproc)
    chainfractions[,]<-sapply(1:nproc,function(j) sapply(1:nproc,function(i) 
        if(i==j){
            1
        }else if(j<i){
            0
        }else{
            Reduce("*",colSums(lam3[1:nproc,1:nproc])[i:(j-1)])
        }
    ))
    
    
    #     chainfractions<-matrix(0,ncol=nproc,nrow=nproc)
    #     for(i in 1:1){
    #         for(j in 1:nproc){
    #             chainfractions[i,j]<-if(j>i){flowfin[i,j]}else{1}*flowfin[j,j]
    #             for(jstar in 1:(i-1)){
    #                 chainfractions[i,j]<-Reduce("*",)
    #             }
    #         }
    #     }
    
    
    
    
    # Share per process which is exported
    lamexport<-colSums(lam3[flows$exportflows[1:nprod],])
    
    # Fate of recycled input to processes
    # The rows give the process into which the input occurs
    # The columns give the process to which the input is distributed
    recfate<-t(t(chainfractions)*lamexport)
    return(recfate=recfate)
}

f_newnutflow<-function(E,S,flows,lambda){
    P<-E$P
    nproc<-S$nproc
    nprod<-S$nprod
    origin<-flows$origin
    target<-flows$target
    goods<-S$goods
    losses<-S$losses
    prows<-S$prows
    sumlosses<-flows$sumlosses
    
    # Fill in input flows
    #sel<-which(origin!=target&target!=0)
    #A<-P[1:nprod,]*lambda
    #for(i in 1:length(sel)){A[sel[i],target[sel[i]]]<--A[sel[i],origin[sel[i]]]}
    A<-P[1:nprod,]
    #Add to A the losses according to flow-allocation
    lbyflow<-allocationbyflow(E,S,flows)
    ladd<-t(t(lbyflow)*sumlosses)
    #Substract now the losses according to current allocation
    lsub<-t(t(lambda)*sumlosses)
    A<-A+ladd-lsub
    #Adjust recycling flows
    for(i in 1:length(sel)){A[sel[i],target[sel[i]]]<--A[sel[i],origin[sel[i]]]}
    
    sel<-which(rownames(l2)%in%goods)
    #print(apply(as.matrix(Aall[sel,][origin[sel]==3,]),2,sum))
    #print(P[sel,][origin[sel]==3,,drop=FALSE])
    V<-matrix(rep(0,nproc**2),ncol=nproc,nrow=nproc)
    V<-t(sapply(1:nproc,function(x) apply(l2[sel,][origin[sel]==x,,drop=FALSE],2,sum)))
    e<-as.vector(as.matrix(-P[prows%in%losses,]))
    rintensity<-t(r)%*%ginv(V)
    eintensity<-t(e)%*%ginv(V)
    
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
            #QQ Is flows really needed to be done or can it be omitted?
            E$P<-systemseparation(E,S,flows,lambda)
            testnut<-nutflowanalysi(E,S,flows)
            
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
            
            #print("chainfractions")
            #print(chainfractions)
            recfateold<-f_recfate(S,flowrec,chainfractions)
            recfate<-f_recfatelambda(E,S,flows,lambda,flowfin)
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
