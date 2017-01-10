library("Matrix")                       # Matrix calculations
library("MASS")                         # Inverted matrix
#Example Heijungs and Suh, 2002, page 16ff

#A: Technology matrix. economic flows: fuel (intermediate) and electricity (demanded)
# Economic flows are flows which come from or go to another process (Definition 1, p. 20)
A<-matrix(c(-2,10,100,0),nrow=2,ncol=2)

# a vector of economic flows
# The economic vector is partitioned into two sets, goods G and wastes W (p. 61)
# a<-c(g,w)

# Muti-functionality
# 3.90 Co-production: Two or more outputs of goods (gi < 0 and gk < 0)
# 3.91 Co-waste treatment: To or more inputs of wastes (wi < 0  and wk< 0)
# 3.92 Recycling: One or more outputs of goods and one or more inputs of wastes

# Definition Recycling (p. 67)
# In LCA, recycling is defined as the situation in which a unit process
# trans-forms a negatively valued product or material (i.e. a waste)
# into a positively valued product or material.

# Definition of Close-Loop Recycling (p. 68)
# We speak of closed-loop recycling when secondary material produced by 
# a recycling process is completely fed back into one of the unit processes 
# of the same product system. When the material is transferred to another 
# product system, we speak of open-loop recycling. 


#B: Intervention matrix. environmental flows (emissions of CO2 and sulfur, use of the resource crude oil)
# Environmental flows come from or go to the environment (Definition 1, p. 20)
B<-matrix(c(1,0.1,0,10,2,-50),nrow=3)
#f: demand vector for economic flows
#   there is no demand for fuel (required to generate elecricity in second process)
#   there is demand of 1000 kWh of electricity
f<-c(0,1000)

#P: process matrix 
# Definition 2 A process matrix P is a set of process vectors, juxtaposed
#to one another. It may be partitioned into a technology matrix A that
#represents the exchanges between processes, and an intervention matrix B
#that represents the exchanges with the environment. (Page 21)
P<-rbind(A,B)

#s scaling vector
s<-ginv(A)%*%f

#g Inventory vector. Environmental flows
#Definition 5 An inventory vector g is a vector of environmental flows.
#The coefficients of this vector represent the amount of these items that a
#system under consideration absorbs or produces. (page 22)
g<-B%*%s
g<-B%*%ginv(A)%*%f

#Lambda: Intensity matrix
Lambda<-B%*%ginv(A)

## The substitution approach (p.41)
## The partitioning methods (p.46)
# - allocation introduced 'manually' into the formulae- this has not been formalized

## The surplus method (p.49)
# - basically considering coproducts as residuals

## Regression approach (section 3.2.5 and 3.5.2 p. 69)
## =========================
A<-matrix(c(-2,10,-1,100,0,50),nrow=3,ncol=2)
f<-c(0,1000,0)
s<-ginv(t(A)%*%A)%*%t(A)%*%f #see formula p. 69, derivation p. 53
# Under certain conditions (see Harville (1997, p.495) and Albert (1972)), 
# the matrix A+ is identical to a matrix that is known as the Moore-Penrose 
# inverse or the pseudoinverse. 
Aplus<-ginv(t(A)%*%A)%*%t(A)
ftilde<-A%*%s

#Discrepancy vector d
d<-ftilde-f

