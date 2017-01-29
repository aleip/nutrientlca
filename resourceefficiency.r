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
    f<-E$f
    nproc<-ncol(P)
    products<-S$products
    nprod<-sum(rownames(P)%in%products)
    recycled<-rownames(P)%in%S$recycled*1
    products<-S$products
    goods<-S$goods
    resources<-S$resources
    residues<-S$residues
    losses<-S$losses
    
    #flagging flows from the process they are generated (0 for imports)
    origin<-sapply(1:nprod,function(x) if(length(which(P[x,]>0))){which(P[x,]>0)}else if(length(which(P[x,]<0))){which(P[x,]<0)}else{0})
    target<-sapply(1:nprod,function(x) if(length(which(P[x,]<0))){mean(which(P[x,]<0))}else{0})
    target[which(target==0&recycled[1:nprod]==1)]<-origin[which(target==0&recycled[1:nprod]==1)]
    
    # Definition chainflows: Flows which are completely produced and absorbed 
    #   in the supply chain whereby the consumption process is at 'higher level'
    #   than the production process
    #chainflows<-sapply(1:nprod,function(x) dochainflows(x,P))
    chainflows<-(target>origin)
    
    # Definition recyclingflows: Flows which are completely produced and absorbed
    #   in the supply chain whereby the consumption process is at 'lower level'
    #   than the production process
    recyflows<-(target<origin)&(target!=0)
    
    # Definition residualflows: Flows which are completely produced and absorbed
    #   in the supply chain whereby production and consumption occurs in the 
    #   same process
    resiflows<-recyflows
    resiflows[]<-FALSE
    resiflows[which(rownames(P)%in%residues)]<-TRUE
    target[resiflows]<-origin[resiflows]
    
    # Defintion mainflows: a demand is quantified
    mainflows<-c(f>0,rep(FALSE,nprod-length(f)))
    
    # Definition exportflows: Flows which leave the current supply chain (co-products)
    expoflows<-origin!=0 & target==0 & !mainflows
    #exportflows<-target==0
    
    
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
    
    return(list(origin=origin,
                target=target,
                
                chainflows=chainflows,
                recyflows=recyflows,
                resiflows=resiflows,
                mainflows=mainflows,
                expoflows=expoflows,
                
                suminputs=suminputs,
                sumgoods=sumgoods,
                sumlosses=sumlosses,
                sumresources=sumresources,
                
                nue=nue)
    )
}

systemseparation<-function(E,S,flows,lambda){
    P<-E$P
    goods<-S$goods
    products<-S$products
    resources<-S$resources
    losses<-S$losses
    nprod<-S$nprod
    nproc<-S$nproc
    f<-E$f
    
    origin<-flows$origin
    target<-flows$target
    sumlosses<-flows$sumlosses
    sumresources<-flows$sumresources
    mainflows<-flows$mainflows
    expoflows<-flows$expoflows
    chainflows<-flows$chainflows
    recyflows<-flows$recyflows
    resiflows<-flows$resiflows
    
    
    # Distribute burden over the goods
    em<-t(t(lambda)*sumlosses)
    emmain<-em*(mainflows|chainflows)
    emcopr<-em*(expoflows)
    
    prmain<-P[1:nprod,]*(mainflows|chainflows)
    prmain[prmain<0]<-0
    prcopr<-P[1:nprod,]*(expoflows)
    #Do not consider recycling flows so far
    inmain<-(prmain+emmain)
    incopr<-(prcopr+emcopr)
    inmain[inmain<0]<-0
    incopr[incopr<0]<-0
    
    # Allocation factors based on input to processes
    inallmain<-colSums(inmain)/colSums(inmain+incopr)
    inallcopr<-colSums(incopr)/colSums(inmain+incopr)
    # Allocate resources based on input shares
    resmain<-sumresources*inallmain
    rescopr<-sumresources*inallcopr
    
    #Residues flows are 
    resimain<-t(t(P[resiflows,])*inallmain)
    resicopr<-t(t(P[resiflows,])*inallcopr)
    
    # Distribute recylcing flows between main and co-processes
    # Problem: if recycling flows are considered as co-product and not as residuals
    #          then system separation and application of the matrix-inversion 
    #          approach leads potentially to problems.
    # Therefore recycling flows HAVE to be considered as residuals and are
    # split at production side acc to share of products
    # and at input side acc to share of input (resources)
    # (there are many alternatives solutions that could also be OK, 
    # important is that the split can reproduce MFA for flow-allocation
    # e.g. recycling flows could be split acc to allocation of burden at production side
    # )
    # ----> maybe this doesn't matter as the inversion is done only for the first rows?
    recmaintarget<-t(t(P[1:nprod,]*recyflows)*inallmain)
    reccoprtarget<-t(t(P[1:nprod,]*recyflows)*inallcopr)
    recmaintarget[recmaintarget>0]<-0
    reccoprtarget[reccoprtarget>0]<-0
    
    
    # Allocation shares based on production
    # (necessarily these allocation shares will be equal to resource shares)
    prallmain<-colSums(prmain)/colSums(prmain+prcopr)
    prallcopr<-colSums(prcopr)/colSums(prmain+prcopr)
    recmainorigin<-t(t(P[1:nprod,]*recyflows)*prallmain)
    reccoprorigin<-t(t(P[1:nprod,]*recyflows)*prallcopr)
    emmainorigin<-
        recmainorigin[recmainorigin<0]<-0
    reccoprorigin[reccoprorigin<0]<-0
    
    recmain<-recmaintarget+recmainorigin
    reccopr<-reccoprtarget+reccoprorigin
    #Distribute remaining emissions [if there was any attached to recycling flows]
    emmain[recyflows]<-t(t(em[recyflows,])*prallmain)
    emcopr[recyflows]<-t(t(em[recyflows,])*prallcopr)
    
    # Distribute chainflows
    prmain[1:nproc,1:nproc]<-sapply(1:nproc,function(j) sapply(1:nproc,function(i) 
        if(target[i]==j) {prmain[i,j]-prmain[i,origin[i]]*inallmain[j]}else{prmain[i,j]}))
    prcopr[1:nproc,1:nproc]<-sapply(1:nproc,function(j) sapply(1:nproc,function(i) 
        if(target[i]==j) {prcopr[i,j]-prmain[i,origin[i]]*inallcopr[j]}else{prcopr[i,j]}))
    
    # Re-constructing the Process matrices
    prmain[recyflows,]<-recmain[recyflows,]
    prcopr[recyflows,]<-reccopr[recyflows,]
    
    if(length(resources)==1){
        prmain<-rbind(prmain,resmain)
        rownames(prmain)[which(rownames(prmain)=="resmain")]<-resources
        prcopr<-rbind(prcopr,rescopr)
        rownames(prcopr)[which(rownames(prcopr)=="rescopr")]<-resources
    }else{stop("More than one resource row - please adapt script")}
    if(length(losses)==1){
        prmain<-rbind(prmain,colSums(emmain))
        rownames(prmain)[nrow(prmain)]<-losses
        prcopr<-rbind(prcopr,colSums(emcopr))
        rownames(prcopr)[nrow(prcopr)]<-losses
        rownames(prcopr)[which(rownames(prcopr)=="rescopr")]<-resources
    }else{stop("More than one resource row - please adapt script")}
    
    # Final step - combine the two matrices
    colnames(prcopr)<-paste(colnames(prmain),"Copr")
    prnew<-cbind(prmain[,1],prcopr[,1])
    prnewn<-c(colnames(prmain)[1],colnames(prcopr)[1])
    for(i in 2:nproc){
        prnew<-cbind(prnew,prmain[,i],prcopr[,i])
        prnewn<-c(prnewn,colnames(prmain)[i],colnames(prcopr)[i])
    }
    colnames(prnew)<-prnewn
    
    newcolnams<-colnames(prnew)[!colnames(prnew)%in%rownames(prnew)]
    #print(newcolnams)
    #print(prnew[,newcolnams])
    if(sum(prnew[,newcolnams])!=0){
        stop(paste("Problem in coproduction processes ",newcolnams,"!\n",
                   "Hint: The processes must have the names of the main products that they produce.\n",
                   "All processes generating co-products are names with including 'Copr' to the name.\n",
                   "For multiple co-products script needs to be adapted..."))
    }else{
        #print(newcolnams)
        prnew<-prnew[,colnames(prnew)[!colnames(prnew)%in%newcolnams]]
    }
    
    
    # Splitting recycling flows from two processes to two processes 
    # gives four individual flows
    #   addrecyflows<-function(matr,recyflows){
    
    matr<-prnew
    addi<-0
    ri<-9
    for(ri in which(recyflows)){
        #    for(i in 9){
        addrows<-prnew[ri,]
        addmatr<-matrix(rep(addrows,4),ncol=length(addrows),nrow=4,byrow=TRUE)
        colnames(addmatr)<-colnames(prnew)
        rownames(addmatr)<-rep(rownames(prnew)[ri],4)
        shari<-addrows[addrows>0]/sum(addrows[addrows>0])
        sharo<-addrows[addrows<0] 
        addmatr[,]<-0
        for(j in (1:length(shari))){
            for(i in (1:length(sharo))){
                puti<-which(colnames(addmatr)==names(sharo)[i])
                puto<-which(colnames(addmatr)==names(shari)[j])
                addmatr[length(shari)*(j-1)+i,puti]<-sharo[i]*shari[j]
                addmatr[length(shari)*(j-1)+i,puto]<--sharo[i]*shari[j]
            }
        }
        matr<-rbind(matr[1:(ri-1+addi),],addmatr,matr[(ri+addi+1):nrow(matr),,drop=FALSE])
        addi<-addi+3
    }

    # Splitting chain flows from one process to two processes 
    # gives two individual flows
    #Not sure yet if this is needed
    #addi<-0
    #ri<-1
    #for(ri in which(chainflows)){
    #    #    for(i in 9){
    #    addrows<-prnew[ri,]
    #    addmatr<-matrix(rep(addrows,2),ncol=length(addrows),nrow=2,byrow=TRUE)
    #    colnames(addmatr)<-colnames(prnew)
    #    rownames(addmatr)<-rep(rownames(prnew)[ri],2)
    #    shari<-addrows[addrows>0]/sum(addrows[addrows>0])
    #    sharo<-addrows[addrows<0] 
    #    addmatr[,]<-0
    #    for(j in (1:length(shari))){
    #        for(i in (1:length(sharo))){
    #            puti<-which(colnames(addmatr)==names(sharo)[i])
    #            puto<-which(colnames(addmatr)==names(shari)[j])
    #            addmatr[length(shari)*(j-1)+i,puti]<-sharo[i]*shari[j]
    #            addmatr[length(shari)*(j-1)+i,puto]<--sharo[i]*shari[j]
    #        }
    #    }
    #    if((ri-1+addi)>0){
    #        print(ri)
    #        print(addi)
    #        print(matr[1:(ri-1+addi),])
    #        matr<-rbind(matr[1:(ri-1+addi),],addmatr,matr[(ri+addi+1):nrow(matr),,drop=FALSE])
    #    }else{
    #        matr<-rbind(addmatr,matr[(ri+addi+1):nrow(matr),,drop=FALSE])
    #    }
    #    #print(matr)
    #    addi<-addi+1
    #}
    
    
        return(P=matr)
}

nutflowanalysi<-function(E,S,flows){
    P<-E$P
    #In case system separation has been done
    S$nproc<-ncol(P)
    nproc<-S$nproc
    
    nprod<-sum(rownames(P)%in%S$products)
    prows<-rownames(P)
    origin<-flows$origin
    resources<-S$resources
    losses<-S$losses
    waste<-S$waste
    goods<-S$goods
    
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
        save(etype,S,file="temp")
        E<-supplyvalues(etype,S)
        
        s<-scaling(E$A,E$B,E$f)
        E$P<-s$P
        E$P<-fillrecycling(E,S)
        
        flows<-flowanalysis(E,S)
        
        # B. Calculation of material flow analysis
        mfa<-nutflowanalysi(E,S,flows)
        
        # C. Calculation of efficiency acc to allocation
        if(dolambda=="byflow"){
            lambda<-allocationbyflow(E,S,flows) 
        }else if(dolambda=="byvalue1"){
            lambda<-allocationbyvalue1(E$P,S$nprod,S$nproc) 
        }else if(dolambda=="byvalue2"){
            lambda<-allocationbyvalue2(E$P,S$nprod,S$nproc) 
        }
        
        E$P<-systemseparation(E,S,flows,lambda)
        # Redo flow analysis for new Process matrix
        flows<-flowanalysis(E,S)
        lca<-nutflowanalysi(E,S,flows)
        
        curex<-list(stype=stype,etype=etype,dolambda=dolambda,lambda=lambda,
                    S=S,E=E,flows=flows,mfa=mfa,lca=lca)
        if(i==1){examp<-curex}else if(i==2){examp<-list(examp,curex)}else{
            examp[[i]]<-curex
        }
    }
    return(examp)
}

comparison<-function(a,what2,which2=NULL,showmfa=FALSE){
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
            mainf<-b$S$products[b$flows$mainflows]
            chainf<-b$S$products[b$flows$chainflows]
            showf<-which(colnames(b$E$P)%in%c(mainf,chainf))
            if(showmfa){
                approach<-"MFA"
                cat(i,approach,round(as.vector(1/b$mfa$rintensity),2),
                    "-",b$stype,"-",b$etype,"\n")
            }
            #else{
            approach<-"LCA"
            b$burden$nueproducts[b$burden$nueproducts==0]<-NA
            cat(i,approach,round(as.vector(1/b$lca$rintensity[showf]),2),
                "-",b$stype,"-",b$etype,"-",b$dolambda,"\n")
            #}
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
