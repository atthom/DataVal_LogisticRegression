library(lattice) 
library(rgl) 
library(akima)
library(plot3D)
library(rgdal)
library(pixmap)
library(gplots)
library(pheatmap)
n = 1000
m0 = 0.5
m1 = 1.1
sd = 1.2

#1.
x1 = rnorm(n,m0,sd)

#2.
x2 = rnorm(n,m1,sd)

#3.
b0 = 0.5
b1 = 1.5

logistic<-function(z){
  return(1.0/(1.0+exp(-z)))
}

generateBinaryResponse<-function(x,b0,b1){
  y = rbinom(n,1,logistic(b0+b1*x))
  return(y)
}

j_funct <- function(xi, yi, m, bc) {
  all_sum =  c(0,0)
  for(i in 1:m) {
    one = (get_pc(bc, xi[i]) - yi[i])* xi[i]
    all_sum = c(all_sum, one)
  }
  return((1/m)*sum(all_sum))
}

gradient<-function(xi,yi,bc){
  
  correction = yi - get_pc(bc, xi)
  print( get_pc(bc, xi))
  
  gradJ = rowSums(xi*correction)
  
  return(gradJ)
}


y1 = generateBinaryResponse(x1,b0,b1)

#4.
y2 = generateBinaryResponse(x2,b0,b1)

#5.
x = c(x1,x2)
xi = t(cbind(rep(1,length(x)),x))
yi = c(y1,y2)

bc0 = 0.3
bc1 = 1.7
bc = c(bc0,bc1)

gamma = 0.0001

M = 200
ccc = c(0,0)

compute_j<-function(yi, beta0, beta1,x) {
  # function that compute J of beta
  return(sum(yi%*%(beta0+beta1*x))-sum(log(1+exp(beta0+beta1*x))))
}

bc_0_all = bc[1]
bc_1_all = bc[2]

all_gradJ = c()

for(t in 1:M){
  #We compute the gradient
  pc = logistic(bc[1]+bc[2]*x)
  gradJ = rowSums(xi%*%(yi-pc))
  cat("gradJ = ",gradJ,"\n")
  all_gradJ = c(all_gradJ, norm(as.matrix(gradJ)))
  

  #We update the estimates 
  bc = bc + gamma*gradJ
  
  #We compute the maximum likelihood function
  J = compute_j(yi, bc[1], bc[2], x)
  cat("J = ",J,"\n")
  ccc = c(ccc, J)
  
  
}

plot(all_gradJ, title(main = "Norme du gradient en fonction de t", xlab = "iterations", ylab = "Norme du gradient"))


bc_0_all = seq(0, 2, 0.1)
bc_1_all = seq(0, 2, 0.1)
ll = length(bc_1_all)
mat = matrix(rep(0, (ll)*(ll)), ncol=ll, nrow=ll)

for(i in 1:(length(bc_1_all))) {
  for(j in 1:(length(bc_1_all))) {
    mat[i,j] =  compute_j(yi, bc_0_all[i], bc_1_all[j], x)
  }
}

colnames(mat) = paste("Beta 0 : ", bc_0_all , sep = "")
rownames(mat) = paste("Beta 1 : ", bc_1_all, sep = "")
pheatmap(mat, cluster_row = FALSE, cluster_col = FALSE, color= heat.colors(254, alpha = 1),
         main = "Heatmap de la fonction de coup (J) en fonction de beta 0 et de beta 1")

cat("estimated beta = ",bc,"\n")

#7.
Prm0 = sum(generateBinaryResponse(x1,bc[1],bc[2]))/n
Prm1 = sum(generateBinaryResponse(x2,bc[1],bc[2]))/n

