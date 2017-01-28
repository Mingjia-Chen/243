getwd()
setwd("~/Google Drive/STAT243/final")
data <- read.csv("wines.csv",sep=",")
data1 <- data[,2:54]
sets <- list(2:7,8:13,14:19,20:24,25:30,31:35,36:39,40:45,46:50,51:54)

cut_data <- function(data,sets){
  for(i in 1:length(sets)){
    data_sub <- data[,sets[[i]]]
    assign(paste("data",i,sep = ""),data_sub,envir = .GlobalEnv)
  }
}


cut_data(data,sets)

##pre_processed
pre_processed <- function(data){
  data_mean <- apply(data,2,mean)
  data_center <- sweep(data,2,STATS = data_mean,FUN = "-")
  data_sum <- apply(data_center,2,function(x) {sqrt(sum(x^2))})
  data_norm <- sweep(data_center,2,STATS = data_sum,FUN = "/")
  data_norm
}
##svd
sum <- 0
a <- NULL
Q <- NULL
for(i in 1:10){
  #name <- sprintf("data%d",3,sep = "")
  l <- length(sets[[i]])
  SVD <- svd(pre_processed(data[,sets[[i]]]))
  #V <- SVD$v
  gamma <- SVD$d[1]
  alpha <- 1/gamma^2
  a[(sum+1):(sum+l)] <- rep(alpha,l)
  #A <- diag(a[(sum+1):(sum+l)])
  #Q <- rbind(Q,solve(A)*V)
  sum <- sum+l
}

A <- diag(a)
M <- diag(rep(1/nrow(data),nrow(data)))
#M <- diag(rep(0.08,12))

#P <- chol(solve(M))
#Q <- chol(solve(A))
#Delta <- solve(P)%*%data%*%solve(t(Q))

#assign(paste("data",i,"norm",sep = ""),svd(pre_processed(data[,sets[[i]]]))$d[1],envir = .GlobalEnv)
#pre <- pre_processed(data1)
#pre <- scale(data1)
X <- pre_processed(data1)
pre <- as.matrix(X)
SVD <- svd(pre)
D <- SVD$d
U <- -SVD$u
V <- SVD$v

alpha_1 <- 1/(D[1])^2
G <- U*D
#A <- diag(rep(0.241,6))
A1 <- sqrt(A)
Q <- solve(A1)%*%V


Q1 = diag(1/sqrt(diag(A)), nrow(A)) %*% V


M1 <- sqrt(M)
P <- solve(M1)%*%U
#t(Q)%*%
P_inv <- t(P)%*%M
Qt_inv <- A%*%Q
Delta <- P_inv %*% pre %*% Qt_inv
F_matrix <- P%*%Delta
#as.matrix(pre)%*%A%*%Q
F_matrix



a <- c(rep(0.241,6),rep(0.239,6),rep(0.275,6),rep(0.273,5),rep(0.307,6),rep(0.302,5),rep(0.417,4),rep(0.272,6),rep(0.264,5),rep(0.309,4))




X <- pre_processed(data1)
X <- as.matrix(X)
s <- X%*%A%*%t(X)
svd <- svd(s)
u <- svd$u
v <- svd$v
d <- svd$d
p <- sqrt(solve(M))%*%u
delta2 <- solve(p)%*%s%*%solve(t(p))
p <- p[,1:11]
delta <- diag(sqrt(diag(delta2)[1:11]))
#q <- sqrt(solve(A))%*%v
q <- t(X)%*%M%*%p%*%solve(delta)

X1 <- p%*%delta%*%t(q)
X1[1:5,1]
X[1:5,1]

##common factor score
f <- p%*%delta




##partial factor score
K <- length(sets)
sum <- 0
for(i in 1:K){
  l <- length(sets[[i]])
  assign(paste("F_",i,sep = ""), K*a[(sum+1)]*X[,(sum+1):(sum+l)]%*%q[(sum+1):(sum+l),])
  sum <- sum+l
}


sum <- 0
dim(q)
ctr_table <- NULL
ctr_obs <- NULL
ctr_var <- NULL
for(i in 1:K){
  l <- length(sets[[i]])
  ctr_var <- rbind(ctr_var,a[(sum+1)]*q[(sum+1):(sum+l),]^2)
  ctr_table <- rbind(ctr_table,apply(a[(sum+1)]*q[(sum+1):(sum+l),]^2,2,sum))
  sum <- sum+l
}

L <- 11

for(i in 1:L){
  
}

lambda <- apply(1/nrow(data1)*f^2,2,sum)
ctr_obs <- sweep(1/nrow(data1)*f^2,2,STATS = lambda, FUN = "/")
