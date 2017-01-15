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
    
    #flagging recycled flows
    recycled<-rep(0,length(products))
    recycled[6:7]<-1
    
    waste<-c()
    goods<-products[!products%in% waste]
    prows<-c(products,resources,losses)
    
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
