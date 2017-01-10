chain1[1,1]<-chain1[1,1]/stages1[1,nflows]
stages1[1,]<-stages1[1,]/stages1[1,nflows]

inputfactors1<-finputfactors(nstages,stages1)
flows1<-fproducts(inputfactors1,product)*stages1
netout1<-apply(fflows(nstages,stages1,product)[,3:nflows],1,sum)
totinputs1<-inputfactors1*product
totoutputs1<-apply(flows1[,2:nflows],1,sum)
nue1<-netout1/totinputs1
extinputs1<-fextinputs(nstages,chain1,totinputs1,product)
stages1a<-stages1*fproducts(infactors1,product)
chain1a<-chain1*product
mylca1<-lcanue(nstages,product,chain1,netout1,extinputs1)

cat("Total input=",totinputs1,"\nExternal input=",extinputs1,
    "\nExport=",stages1a[,3],"\nLosses=",stages1a[,2],
    "\nResidues=",stages1a[,5],"\nRecyled=",stages1a[,4],"\nProducts=",products(infactors1,product),
    "\nNUE for each stage=",nue1,"\nLife-cycle NUE=",mylca1,
    "\nTotal inputs=",totinputs1," Total outputs=",totoutputs1,
    "\n\n")
