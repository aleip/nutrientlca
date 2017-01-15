upperTriangle<-function(x){
    n<-min(ncol(x),nrow(x))
    d<-sapply(1:n,function(j) sapply(1:n,function(i) if(j<i){0}else{1}))
    m<-max(nrow(x),ncol(x))-n
    f<-matrix(0,ncol=n,nrow=m)
    if(ncol(x)>nrow(x)){
        g<-cbind(d,t(f))
    }else{
        g<-rbind(d,f)
    }
    return(g)
}
