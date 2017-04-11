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
  return(1/(1+exp(-z)))
}

generateBinaryResponse<-function(x,b0,b1){
  y = logistic(b0+b1*x)
  y = floor(y+0.5)
  return(y)
}

y1 = generateBinaryResponse(x1,b0,b1)

#4.
y2 = generateBinaryResponse(x2,b0,b1)
