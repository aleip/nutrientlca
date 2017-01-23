burden$nueproducts[burden$nueproducts==0]<-NA
refburden$nueproducts[refburden$nueproducts==0]<-NA

if(!noPprint){
    cat("Change of the Process matrix:\nOriginal matrix:\n")
    print(refE$P)
    cat("\nChanged matrix (note that the process are scaled to generate the same final demand of processed food):\n")
    print(E$P)
    cat("\n")
}

cat("The Nutrient Use Efficiency of the MFA approach\n")
cat("MFA old",round(as.vector(1/refnutflow$rintensity),2),"\n")
cat("MFA new",round(as.vector(1/nutflow$rintensity),2),"\n")

cat("\nThe Nutrient Use Efficiency of the LCA approach (",
    "allocation ",dolambda,")\n")
cat("LCA old",round(apply(refburden$nueproducts,2,mean,na.rm=TRUE),2),"\n")
cat("LCA new",round(apply(burden$nueproducts,2,mean,na.rm=TRUE),2),"\n")

#print(nutflow$lcanue)
#cat("\nResource Use\n")
#cat("MFA old",round(as.vector(colSums(refburden$finproducts)*refnutflow$rintensity),0),"\n")
#cat("MFA new",round(as.vector(colSums(burden$finproducts)*nutflow$rintensity),0),"\n")
#cat("\nLCA old",round(colSums(refburden$finproducts)/apply(refburden$nueproducts,2,mean,na.rm=TRUE),0),"\n")
#cat("LCA new",round(colSums(burden$finproducts)/apply(burden$nueproducts,2,mean,na.rm=TRUE),0),"\n")
