library("Matrix")                       # Matrix calculations
library("MASS")                         # Inverted matrix

rm(list=objects())

source("genfunctions.r")
source("resourceefficiency.r")

supplychainsimple<-function(type){
    if(type=="default"){
        processes<-c("Feed",
                     "Livestock",
                     "Food")
        products<-c("Feed","Livestock","Food",
                    "Feed Copr","Livestock Copr",
                    "Feed Residues","Feed Residues",
                    "Waste","Livestock Recy","Food Recy")
        resources<-c("N fertilizer")
        losses<-c("N losses")
        interventions<-c(resources,losses)
        nenv<-length(resources)+length(losses)
        nproc<-length(processes)
        nprod<-length(products)
        
        #flagging recycled flows
        recycled<-c("Livestock Recy","Food Recy")
        residues<-c("Feed Residues")
        
        waste<-c("waste")
        goods<-products[!products%in% waste]
        prows<-c(products,resources,losses)
    }
    return(list(processes=processes,
                nproc=nproc,
                products=products,
                nprod=nprod,
                resources=resources,
                residues=residues,
                losses=losses,
                waste=waste,
                goods=goods,
                interventions=interventions,
                nenv=nenv,
                recycled=recycled,
                prows=prows))
}
supplyvalues<-function(type,S){
    nproc<-S$nproc
    nprod<-S$nprod
    nenv<-S$nenv
    products<-S$products
    
    # Technology matrix
    A<-matrix(0,ncol=nproc,nrow=nprod)
    A[,1]<-c(100,  0,   0,40,0,10,-10,0,-10,-2)
    # The vector is constructed such that consumption of recyling flows from 
    # other processes are not explicit to leave cross-dependencies flexible.
    # They are indicated with NA.
    A[,1]<-c( 25,  0,   0,20,0, 5, -5, 0,NA,NA)
    A[,1]<-c( 50,  0,   0,40,0,10,-10, 0,NA,NA)
    A[,2]<-c(-50, 20,   0, 0,5, 0,  0, 0,10, 0)
    A[,3]<-c(  0,-20,12.5, 0,0, 0,  0, 0, 0, 2)
    
    # Intervention matrix
    
    # Required resources are coming from the resources (here mineral fertilizer)
    # NA will be later replaced and the resources adjusted accordingly
    # This does NOT affect the calculated NUE per process or the lifecycle NUE based on material flow analysis.
    B<-matrix(0,ncol=nproc,nrow=nenv)
    B[,1]<-c(-56,11)
    B[,1]<-c(-112,22)
    B[,2]<-c(   0,15)
    B[,3]<-c(   0,5.5)
    
    # Demand
    f<-c(0,0,12.5) #Demand Vector required for main food chain
    
    if(sum("largemanurerecycling"%in%type)){
        A[,2]<-c(-50, 10,   0, 0,5, 0,  0,0,20, 0)
    }
    
    if(sum(type%in%c("nofoodres","norecycling"))){
        temp<-which(products=="Food Recy")
        # Reduce this to 0.01 and add it to final product
        A[3,3]<-A[3,3]+A[temp,3]
        A[temp,3]<-0
    }
    if(sum(type%in%c("manureexport","norecycling"))){
        tempman<-which(products=="Livestock Recy")
        # Manure not recycled by exported (same value)
        # Substitution by mineral fertilizer 'automatic'
        #B[1,1]<-B[1,1]-A[tempman,2]
        A[tempman,2]<-0
    }
    if(sum(type%in%c("morefeedrecycling"))){
        A[,1]<-c( 25,  0,   0,40,0,35,-35, 0,NA,NA)
    }
    if(sum(type%in%c("nofeedrecycling"))){
        A[,1]<-c( 50,  0,   0,40,0,0,0, 10,NA,NA)
        B[,1]<-c(-122,22)
    }
    if(sum(type%in%c("cropresexport"))){
        A[,1]<-c( 50,  0,   0,50,0,0,0, 0,NA,NA)
        B[,1]<-c(-122,22)
    }
    if(sum(type%in%c("food2feed"))){
        A[,1]<-c( 50,  0,   0,40,0,10,-10, 0,NA, 0)
        A[,2]<-c(-48, 20,   0, 0,5, 0,  0, 0,10,-2)
    }
    
    rownames(A)<-S$products
    rownames(B)<-S$interventions
    colnames(A)<-S$processes
    colnames(B)<-S$processes
    
    return(list(A=A,B=B,f=f))
    
}

allocationbyvalue1<-function(P,nprod,nproc){
    #Define relative 'value' (or whatever the basis for allocation is) over the goods produced in each process
    lam2<-matrix(0,ncol=nproc,nrow=nprod)
    # Main products 'double value' per unit of N
    lam2[1,1]<-P[1,1]*1.2
    lam2[2,2]<-P[2,2]*1.2
    lam2[3,3]<-P[3,3]*1.2
    # Exported products 'normal value
    lam2[4,1]<-P[4,1]*1
    lam2[5,2]<-P[5,2]*1
    # Recycled products backwards half value
    lam2[9,2]<-P[9,2]*0.8
    lam2[10,3]<-P[10,3]*0.8
    # Recycling to same process no burden
    lam2[6,1]<-P[6,1]*0
    # Waste leaves burden
    lam2[8,1]<-P[8,1]*0
    
    lam2<-t(t(lam2)/colSums(lam2))
    
    # lam3 is the allocation to the products that are not flowing back into the supply chain
    # (thus exports or connecting flows)
    #lam3<-lam2*(!chainflows)
    #lam3<-t(t(lam3)/colSums(lam3))
    
    return(lambda=lam2)
    
}
allocationbyvalue2<-function(P,nprod,nproc){
    #Define relative 'value' (or whatever the basis for allocation is) over the goods produced in each process
    lam2<-matrix(0,ncol=nproc,nrow=nprod)
    # Main products 'double value' per unit of N
    lam2[1,1]<-P[1,1]*5
    lam2[2,2]<-P[2,2]*5
    lam2[3,3]<-P[3,3]*5
    # Exported products 'normal value
    lam2[4,1]<-P[4,1]*1
    lam2[5,2]<-P[5,2]*1
    # Recycled products backwards half value
    lam2[9,2]<-P[9,2]*0.8
    lam2[10,3]<-P[10,3]*0.8
    # Recycling to same process no burden
    lam2[6,1]<-P[6,1]*0
    # Waste leaves burden
    lam2[8,1]<-P[8,1]*0
    
    lam2<-t(t(lam2)/colSums(lam2))
    
    # lam3 is the allocation to the products that are not flowing back into the supply chain
    # (thus exports or connecting flows)
    #lam3<-lam2*(!chainflows)
    #lam3<-t(t(lam3)/colSums(lam3))
    
    return(lambda=lam2)
    
}

