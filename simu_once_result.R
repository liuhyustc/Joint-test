########## Testing results ##############
source("joint_test.R")
source("pchisqsum2.R")
library(foreach)
library(doParallel)
library(dplyr)

simu1 <- function(i,index.case,index.control,all_variables,fd,fy){
  n1 <- n0 <- 2000
  set.seed(i+1)
  case1 <- sample(index.case,size = n1)
  control1 <- sample(index.control,size = n0)
  
  D <- all_variables$D[c(case1,control1)]
  Y <- all_variables$Y[c(case1,control1)]
  X <- all_variables$X[c(case1,control1),]
  G <- all_variables$G[c(case1,control1),]
  
  try({ res <- GetJoint_minp(D,Y,X,G,fd,fy,K=c(seq(0,0.99999,0.1),0.99))
      minp_pros <- res$pros$pactual
      minp_retr <- res$retr$pactual
      minp_EB <- res$EB$pactual
      minp_Ts1 <- res$Ts1$pactual
      minp_Ts2 <- res$Ts2$pactual})
  return(c(minp_pros,minp_retr,minp_EB,minp_Ts1,minp_Ts2))
}

cl <- makeCluster(64)
registerDoParallel(cl)
# c(log(1.2),log(1.5))
# c("null","alter")

# Rho <- c(0,0.1,0.3,0.5,0.7,0.9)
Rho <- c(0,0.2,0.6,0.9)
M <- 1000
for(delta in c(4) ){
  for(test in c("null","alter")){
    for(alphaY in c(log(1.5) )){
      result <- matrix(0,nrow = 6,ncol = length(Rho))
      for(i in 1:length(Rho)){
        rho <- Rho[i]
        load(paste0( "joint_test/",test,"/alphaY", exp(alphaY),"_rho", rho, "_delta", delta, ".RData"))
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
        # set.seed(1234)
        
        system.time(mydata <- foreach(i=1:M,.combine = cbind,
                                      .packages = c("CompQuadForm","aSPU","MASS"))%dopar%
                      simu1(i,index.case = index.case,index.control = index.control
                            ,all_variables =  all_variables, fd = fd, fy = fy)) %>% print()
        
        
        # print("done")
        result[,i] <- c(rho,round(rowMeans(mydata < 0.05),3))
        rownames(result) <- c("rho","pros","retr","EB","Ts1","Ts2")
      }
      
      write.table(result, file=paste0( "joint_test/result/",test,"_alphaY" ,exp(alphaY),"delta",delta, ".txt"))
    }
  }
}


stopCluster(cl)


