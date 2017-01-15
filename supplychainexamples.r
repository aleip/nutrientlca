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
    
    return(list(A=A,B=B,f=f))
    
}
examplvar1<-function(nproc,nprod,nenv){
    
    # Based on exampleaimable
    # Only difference: last process 
    
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
    
    return(list(A=A,B=B,f=f))
    
}
