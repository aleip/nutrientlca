#R code to calculate the life-cycle-NUE: 
library("Matrix")                       # Matrix calculations
library("MASS")                         # Inverted matrix
u = matrix(c(10, 10, 2, 50, 0,   0, 0, 20, 0), nrow=3, ncol=3, byrow=F)
#Resource
r = matrix (c(100, 0, 0), nrow=3, ncol=1, byrow=F)
#change in stock
s = matrix (c(0, 0, 0), nrow=3, ncol=1, byrow=F) # Matrix -s' is filled with a negative sign as it is an output.
#Losses
w = matrix (c(22, 15, 5.5), nrow=3, ncol=1) 
#Product
v = matrix (c(100,0,0,0,35,0,0,0,14.5 ), nrow=3, ncol=3, byrow=F)
#Check validity
ut=t(u)
verf_1 = rowSums(ut, na.rm=F)+t(r)
verf_2 = rowSums(v, na.rm=F)+t(s)+t(w)
nue = ((rowSums(v, na.rm=F)+t(s))/(rowSums(ut, na.rm=F)+t(r)))*100
#nue = (82%, 70%, 72%)
gs = matrix (c(1:9*0),nrow=3, ncol=3, byrow=T)
rf = t(r) %*% ginv (t(v) - u  + gs ) #C
life_cycle_nue = (1/rf[1,3])*100
#life_cycle_nue = 36%
