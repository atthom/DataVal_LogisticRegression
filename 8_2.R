n = 10
m0 = 0.5
m1 = 1.1
sd = 1.2

#1.
x1 = rnorm(n,m0,sd)

#2.
x2 = rnorm(n,m1,sd)

#3.
b0 = 0.3
b1 = 1.7

logistic<-function(z){
  return(1.0/(1.0+exp(-z)))
}

generateBinaryResponse<-function(x,b0,b1){
  y = rbinom(n,1,logistic(b0+b1*x))
  return(y)
}

y1 = generateBinaryResponse(x1,b0,b1)

#4.
y2 = generateBinaryResponse(x2,b0,b1)

#5.
x = c(x1,x2)
xi = t(cbind(rep(1,length(x)),x))
#print(xi)
yi = c(y1,y2)
#print(yi)
bc0 = 1
bc1 = 1
bc = c(bc0,bc1)

gamma = 0.01

M = 1
for(t in 1:M){
    #We compute the gradient
    #print("bc")  
    #print(bc)
    #print("xi")  
    #print(xi)
    #print(bc[1]+bc[2]*x)
    pc = logistic(bc[1]+bc[2]*x)
    #cat("temp = ",temp,"\n")
    #print("yi")
    #print(yi)
    #print("pc")
    #print(pc)
    #print("temp")
    #print(temp)
    #print("xi*temp")
    #print(xi*temp)
    #print("rowSums(xi*temp)")
    #print(rowSums(xi*temp))
    #print(dim(xi))
    #print(dim(temp))
    gradJ = rowSums(xi*(yi-pc))
    cat("gradJ = ",gradJ,"\n")
    
    #gradJb0 = sum(x[0:t]*(y-pc0))
    #gradJb1 = sum(x[0:t]*(y-pc1))
    #gradJb0 = sum(y*x)-sum((exp(bc0)/(1+exp(bc0)))*x)
    #gradJb1 = sum(y*x)-sum((exp(bc1*x)/(1+exp(bc1*x)))*x)
    #cat("gradJb0 = ",gradJb0,"\n")
    #cat("gradJb1 = ",gradJb1,"\n")
  
    #We update the estimates 
    bc = bc + gamma*gradJ
    #bc0 = bc0 + gamma*gradJb0
    #bc1 = bc1 + gamma*gradJb1
    
    #We compute the maximum likelihood function
    J = sum(yi*bc*xi)-sum(log(1+exp(bc*xi)))
    cat("J = ",J,"\n")
}

print(bc)
