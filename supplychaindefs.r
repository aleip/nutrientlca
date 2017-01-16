library("Matrix")                       # Matrix calculations
library("MASS")                         # Inverted matrix

rm(list=objects())

source("genfunctions.r")
source("resourceefficiency.r")

supplychainsimple<-function(type){
    if(type=="default"){
        processes<-c("Feed production",
                     "Livestock production",
                     "Food production")
        products<-c("feed","anim","food","feedco","animco","cropreso","cropresi","manure","foodres")
        resources<-c("N fertilizer")
        losses<-c("N losses")
        nenv<-length(resources)+length(losses)
        nproc<-length(processes)
        nprod<-length(products)
        
        #flagging recycled flows
        recycled<-rep(0,length(products))
        recycled[6:7]<-1
        
        waste<-c()
        goods<-products[!products%in% waste]
        prows<-c(products,resources,losses)
    }
    return(list(processes=processes,
                nproc=nproc,
                products=products,
                nprod=nprod,
                resources=resources,
                losses=losses,
                waste=waste,
                goods=goods,
                nenv=nenv,
                recycled=recycled,
                prows=prows))
}
example<-function(type,S){
    nproc<-S$nproc
    nprod<-S$nprod
    nenv<-S$nenv
    products<-S$products
    
    # Technology matrix
    A<-matrix(0,ncol=nproc,nrow=nprod)
    A[,1]<-c(100,  0,   0,40,0,10,-10,-10,-2)
    # The vector is constructed such that consumption of recyling flows from 
    # other processes are not explicit to leave cross-dependencies flexible.
    # They are indicated with NA.
    A[,1]<-c( 25,  0,   0,20,0, 5, -5, NA,NA)
    A[,1]<-c( 50,  0,   0,40,0,10,-10, NA,NA)
    A[,2]<-c(-50, 20,   0, 0,5, 0,  0, 10, 0)
    A[,3]<-c(  0,-20,12.5, 0,0, 0,  0,  0, 2)
    
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
        A[,2]<-c(-50, 10,   0, 0,5, 0,  0,20, 0)
    }
    
    if(sum(type%in%c("nofoodres","norecycling"))){
        temp<-which(products=="foodres")
        # Reduce this to 0.01 and add it to final product
        A[3,3]<-A[3,3]+A[temp,3]-0.01
        A[temp,3]<-0.01
    }
    if(sum(type%in%c("manureexport","norecycling"))){
        tempman<-which(products=="manure")
        # Manure not recycled by exported (same value)
        # Substitution by mineral fertilizer 'automatic'
        #B[1,1]<-B[1,1]-A[tempman,2]
        A[tempman,1]<-0
    }
    
    return(list(A=A,B=B,f=f))
    
}

