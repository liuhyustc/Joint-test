####### Estimation results #########
source("joint_test.R")
source("pchisqsum2.R")
library(foreach)
library(doParallel)
library(dplyr)

simu1 <- function(i,index.case,index.control,all_variables,fd,fy){
  n1 <- n0 <- 2000
  case1 <- sample(index.case,size = n1)
  control1 <- sample(index.control,size = n0)
  
  D <- all_variables$D[c(case1,control1)]
  Y <- all_variables$Y[c(case1,control1)]
  X <- all_variables$X[c(case1,control1),]
  G <- all_variables$G[c(case1,control1),]
  est <- GetJoint_estimation(D,Y,X,G,fd,fy)
  
  # return(c(minp_pros,minp_retr,kopt_pros,kopt_retr))
  return(est)
}

cl <- makeCluster(64)
registerDoParallel(cl)
# c(log(1.2),log(1.5))
# c("null","alter")

# Rho <- c(0,0.1,0.3,0.5,0.7,0.9)
M <- 1000
Rho <- c(0,0.2,0.6,0.9)

for(delta in log(c(1,1.2,1.5,1.8))){
    for(alphaY in c(log(1.5) )){
      for(i in 1:length(Rho)){
        rho <- Rho[i]
        result <- matrix(0,nrow = 9,ncol = 7)
        load(paste0( "joint_test/null/alphaY", exp(alphaY),"_rho", rho, "_delta", exp(delta), ".RData"))
        D <- data$D # disease phenotype
        Y <- data$Y # secondary phenotype
        G <- data$G # genotype data
        X <- data$X # covariates and environment variable
        para <- data$para # true parameters except intercept values, list
        fd <- para$fd
        fy <- para$fy
        index.case <- which(D==1)
        index.control <- which(D==0)
        all_variables <- list(D=D,Y=Y,X=X,G=G)
        
         
        
        # set.seed(3407)
        set.seed(1234)
        
        system.time(mydata <- foreach(i=1:M,.combine = cbind,
                                      .packages = c("CompQuadForm","aSPU","MASS"))%dopar%
                      simu1(i,index.case = index.case,index.control = index.control
                            ,all_variables =  all_variables, fd = fd, fy = fy)) %>% print()
        
        
        # print("done")
        Glm_est <- matrix(unlist(mydata[1,]),ncol = M,nrow=9,byrow = F)
        Glm_sd <- matrix(unlist(mydata[2,]),ncol = M,nrow=9,byrow = F)
        my_est <- matrix(unlist(mydata[3,]),ncol = M,nrow=9,byrow = F)
        my_sd <- matrix(unlist(mydata[4,]),ncol = M,nrow=9,byrow = F)
        
        true_para <- c(para$alpha0,para$alphaX,para$beta0,para$betaX,para$alphaY)
        result[,1] <- true_para
        result[,2] <- rowMeans(Glm_est)-true_para
        result[,3] <- apply(Glm_est, 1, sd)
        result[,4] <- rowMeans(Glm_sd)
        result[,5] <- rowMeans(my_est)-true_para
        result[,6] <- apply(my_est, 1, sd)
        result[,7] <- rowMeans(my_sd)
        colnames(result) <- c("true_para","Bias","SE","SEE","Bias","SE","SEE")
        
        write.table(round(result,3), file=paste0( "2024.1.19/est/alphaY",exp(alphaY),"_rho" ,rho,"_delta",exp(delta), ".txt"))
        
      }
    }
}
stopCluster(cl)
