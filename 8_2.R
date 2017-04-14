n = 1000
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

#On veut generer les Y tel que logit(Pr(Y=1|X=x))=B0+B1*x
#Donc Pr(Y=1|X=x) = logistic(B0+B1*x)
generateBinaryResponse<-function(x,b0,b1){
  y = rbinom(n,1,logistic(b0+b1*x))
  return(y)
}

y1 = generateBinaryResponse(x1,b0,b1)

#4.
y2 = generateBinaryResponse(x2,b0,b1)

#5.
J<-function(x,y,b0,b1){
  return(sum(y%*%(b0+b1*x))-sum(log(1+exp(b0+b1*x))))
}

x = c(x1,x2)
xi = t(cbind(rep(1,length(x)),x))
yi = c(y1,y2)

bc0 = 0.3
bc1 = 1.7
bc = c(bc0,bc1)

gamma = 0.0001

M = 1000
for(t in 1:M){
  #We compute the gradient
  pc = logistic(bc[1]+bc[2]*x)
  gradJ = rowSums(xi%*%(yi-pc))
  cat("gradJ = ",gradJ,"\n")

  #We update the estimates 
  bc = bc + gamma*gradJ
  
  #We compute the maximum likelihood function
  J = sum(yi%*%(bc[1]+bc[2]*x))-sum(log(1+exp(bc[1]+bc[2]*x)))
  cat("J = ",J,"\n")
}

cat("estimated beta = ",bc,"\n")

#7.
Prm0 = sum(generateBinaryResponse(x1,bc[1],bc[2]))/n
Prm1 = sum(generateBinaryResponse(x2,bc[1],bc[2]))/n

#8.
exp(Prm0-Prm1)
