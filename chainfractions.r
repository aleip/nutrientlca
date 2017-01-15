chainfractions<-matrix(0,ncol=nproc,nrow=nproc)
for(j in nproc:1){
    # Loop over all last-iut-one origins for the target
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

finfractions<-matrix(0,ncol=nproc,nrow=nproc)
for(i in 1:nproc){
    
    #if(i==1){finfractions[i,]<-chainfractions[i,]}else{
       # finfractions[i,]<-
    #            # The non-recycled part is distribured as in defined in chainfractions along the chain
    #            flowrec[i,i]*chainfractions[i,]
    #    print(finfractions[i,])
        for(ii in 1:max(1,(i-1))){
            cat(ii,": ",finfractions[i,],"\n")
            cat(ii,": ",flowrec[i,ii],finfractions[ii,],"\n")
            finfractions[i,]<-finfractions[i,]+
                
                # The recycled part enters in a previous step
               flowrec[i,ii]*finfractions[ii,]
            
            cat(i,ii,chainfractions[i,],"\n",
                flowrec[i,i],chainfractions[i,],"\n",
                finfractions[i,],"\n")
        }
        finfractions[i,]<-finfractions[i,]+
                    # The non-recycled part is distribured as in defined in chainfractions along the chain
                    flowrec[i,i]*chainfractions[i,]
        
    #}
    
    
}
#print(flowrec)
#print(chainfractions)
#print(finfractions)