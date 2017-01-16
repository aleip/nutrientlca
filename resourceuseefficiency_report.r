nueproducts<-burden$nueproducts
nueproducts[nueproducts==0]<-NA

cat("Nutrient Use Efficiency\n")
cat("MFA",as.vector(1/nutflow$rintensity),"\n")
cat("LCA",apply(nueproducts,2,mean,na.rm=TRUE),"\n")
#print(nutflow$lcanue)
cat("Resource Use\n")
cat("MFA",as.vector(colSums(burden$finproducts)*nutflow$rintensity),"\n")
cat("LCA",colSums(burden$finproducts)/apply(nueproducts,2,mean,na.rm=TRUE),"\n")
