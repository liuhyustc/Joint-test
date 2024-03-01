args=commandArgs(trailingOnly=TRUE)
N = as.numeric(args[1])  # population
nSNP0 = as.numeric(args[2]) # of casual SNPs
nSNP = as.numeric(args[3])  # of non-casual SNPs
alphaY = log(as.numeric(args[4])) # Y -> D
rho = as.numeric(args[5]) #AR1(rho) between casual SNPs and non-causal SNPs
delta = log(as.numeric(args[6])) # the degree of correlation between E and causal SNPs
test = args[7] ## null or alter
# set = args[6]  ## null--joint/main
# eff = log(as.numeric(args[7]))


# genarate parameters
# N <- 2e6
# nSNP0 <- 5  # generate alphaG and betaG
# nSNP <- 30
# alphaY <- log(1.5)
# rho <- 0
# delta <- log(1)
# test <- "null"

######################## Generate data ##########################################

set.seed(77)
library(Matrix)
library(MASS)
if(test=="null"){
  alphaG <- rep(0, nSNP0)
  betaG <- rep(0, nSNP0)
}else{
  # alphaG <- rnorm(nSNP0,0,0.5) # tau1=1
  # betaG <- rnorm(nSNP0,0,0.5) # tau2=2
  #alphaG <- log(c(1.4,1/1.4,1.2,1/1.2,0.8))
  #alphaG <- rep(0, nSNP0)
  #betaG <- rep(0, nSNP0)
  #betaG <- log(c(2,1/3,3,1/2,0.8))
  alphaG <- log(runif(nSNP0,0.75,1.6))
  # log(runif(5,0.75,1.6))
  betaG <- log(c(1.5,1/1.5,1.2,1/1.2,0.7))
  # alphaG <- log(runif(nSNP0,0.5,2))
  # betaG <- log(c(1.5,1/1.5,1.2,1/1.2,1.8))
}
# alphaX <- c(0.5,0.05,0.1,0.25)
# betaX <- c(0.07,-0.05,0.2,-0.1)
alphaX  <- c(-0.4, -0.13, 0.25) # sex,age,E
betaX <- c(0.2, 0.15, -0.08)
#################### sample size and prevalence and other parameters settup #########
fd <- 0.01
fy <- 0.05
MAF0s <- runif(nSNP0, 0.01, 0.05)
MAFs <- runif(nSNP, 0.01, 0.05)

para <- list(alphaX=alphaX,alphaG=alphaG,alphaY=alphaY,betaX=betaX,betaG=betaG,
             N=N,nSNP0=nSNP0, nSNP=nSNP,MAF0s = MAF0s, MAFs = MAFs, fd=fd,fy=fy,rho=rho,delta=delta)

# para <- list(alphaX=alphaX,alphaG=alphaG,alphaY=alphaY,betaX=betaX,betaG=betaG,
#              N=N,nSNP0=nSNP0, nSNP=nSNP,MAF0s = MAF0s, MAFs = MAFs, fd=fd,fy=fy,rho=rho,eta=eta)

####################### generate population ############
simRareSnp <- function(para){
  alphaX <- para$alphaX
  alphaG <- para$alphaG
  alphaY <- para$alphaY
  betaX <- para$betaX
  betaG <- para$betaG
  N <- para$N
  nSNP0 <- length(alphaG)
  nSNP <- para$nSNP
  fd <- para$fd
  fy <- para$fy
  MAF0s <- para$MAF0s
  MAFs <- para$MAFs
  rho <- para$rho
  delta <- para$delta
  
  
  #W0 <- matrix(diag(dbeta(MAF0s, 1, 25)), ncol = nSNP0) # wu's weights
  # W0 <- matrix(diag(1/sqrt(MAF0s * (1-MAF0s))),ncol = nSNP0)
  R0 <- matrix(1, nrow=nSNP0, ncol=nSNP0)
  for(i in 1:nSNP0)
    for(j in 1:nSNP0)
      if (i!=j) R0[i, j] <- rho^(abs(i-j))
  
  R0noise <- matrix(1, nrow=nSNP, ncol=nSNP)
  for(i in 1:nSNP)
    for(j in 1:nSNP)
      if (i!=j) R0noise[i, j] <- rho^(abs(i-j))

  D <- Y <- numeric(length = N)
  ## generate covariates
  X1 <- rbinom(N, size = 1, prob = 0.5) # sex
  X2 <- rnorm(N ,mean = 30, sd = 9.4) - 30 # age 
  #X3 <- rnorm(N, mean = 169, sd = 9.4)  -169 # height
  ## generate causal SNPs
  eigen.R <- eigen(R0, symmetric = T)
  R1 <- eigen.R$vectors%*% diag(sqrt(eigen.R$values))
  cutoff0 <- qnorm(MAF0s)
  
  G0 <- matrix(NA, nrow = N, ncol = nSNP0)
  Env <- rep(NA, N)   # environment factor
  for (i in 1:N) {
    # all causal SNPs
    r0 <- rnorm(nSNP0, 0, 1)
    r1 <- R1%*%r0 
    r2 <- ifelse(r1<cutoff0, 1, 0)
    
    r0 <- rnorm(nSNP0, 0, 1)
    r1 <- R1%*%r0 
    r3 <- ifelse(r1<cutoff0, 1, 0)
    
    r4 <- r2+r3 ## RV
    G0[i, ] <- r4
    # Env[i] <- sum((r4-MAF0s*2)^2)*delta + rnorm(1, 0, 1)
    Env[i] <- sum(r4)*delta + rnorm(1, 0, 1)
  }
  ## generate neutral SNPs:
  if(nSNP > 0){
    eigen.Rnoise <- eigen(R0noise, symmetric = T)
    R1noise <- eigen.Rnoise$vectors%*%diag(sqrt(eigen.Rnoise$values)) # sqrt(R0noise)
    cutoff <- qnorm(MAFs)
    
    Gnoise <- matrix(0, nrow=N, ncol = nSNP)
    for(i in 1:N){
      r0 <- rnorm(nSNP, 0, 1) #: X0 ~ MVN(0, I)
      r1 <- R1noise %*% r0   #: X1 ~ MVN(0, R)
      r2 <- ifelse(r1 < cutoff, 1, 0)
      
      r0 <- rnorm(nSNP, 0, 1) #: X0 ~ MVN(0, I)
      r1 <- R1noise %*% r0   #: X1 ~ MVN(0, R)
      r3 <- ifelse(r1 < cutoff, 1, 0)
      
      r4 <- r2+ r3
      Gnoise[i, ] <- r4
    }
    Gall <- cbind(G0, Gnoise)
  }else Gall <- G0
  # Env <- (rowSums(G0)-2*MAF0s)*delta + (rowSums(Gnoise)-2*MAFs)*delta+ rnorm(N,0,0.5)
  # generate phenotypes 
  
  X <- cbind(X1,X2,Env)
  beta0 <- Getbeta0(X, G0, para = para, fy = fy)
  PrY <- plogis(rep(beta0,N)+cbind(X,G0)%*%(c(betaX,betaG) ))
  Y <- rbinom(N,1,PrY)
  alpha0 <- Getalpha0(X,G0,Y, para = para, fd = fd)
  PrD <- plogis(rep(alpha0,N)+ cbind(X,G0,Y)%*%(c(alphaX,alphaG,alphaY) ))
  mean(plogis(rep(alpha0,N)+ cbind(X,G0,Y)%*%(c(alphaX,alphaG,alphaY) ))[which(Y==1)] )
  D <- rbinom(N,1,PrD) # sample D according to PrD
  length(intersect(which(D==1),which(Y==1)))
  
  para <- list(alpha0 = alpha0, alphaX = alphaX, alphaG =alphaG, alphaY = alphaY, beta0 = beta0
               ,betaX = betaX, betaG = betaG, nSNP0=nSNP0, nSNP=nSNP,N=N,fd=fd,fy=fy,delta=delta)

  data <- list(D=D, Y=Y, G=Gall, X=cbind(X1,X2,Env),para=para)
  save(data, file=paste0( "joint_test/",test,"/alphaY", exp(alphaY),"_rho", rho, "_delta", exp(delta), ".RData"))
 
}

#############compute beta0 according to fy###############
# Covs <- list(X = cbind(X1,X2,X3,Env), G = Gall , Y =Y)
# alpha <- list(alphaX = alphaX, alphaG = alphaG, alphaY = alphaY)
# beta <- list(betaX = betaX, betaG = betaG)
Getbeta0  <- function(X, G0, para, fy){
  betaX <- para$betaX
  betaG <- para$betaG
  
  zb <- as.vector(cbind(X,G0)%*%(c(betaX,betaG)))
  fn <- function(beta0){
    p <- plogis(beta0 + zb)
    return(mean(p)-fy)
  }
  beta0 <- uniroot(fn, c(-100, 100))$root
  return(beta0)
}

Getalpha0  <- function(X, G0, Y, para, fd){
  alphaX <- para$alphaX
  alphaG <- para$alphaG
  alphaY <- para$alphaY
  
  zb <- as.vector(cbind(X,G0,Y)%*%(c(alphaX,alphaG,alphaY)))
  fn <- function(alpha0){
    p <- plogis(alpha0 + zb)
    return(mean(p)-fd)
  }
  alpha0 <- uniroot(fn, c(-100, 100))$root
  return(alpha0)
}


simRareSnp(para)
rm(list = ls())





  





