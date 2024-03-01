source("pchisqsum2.R")
source("Ind_test.R")
################# prospective likelihood method ##########

################check population data : meet the design########
# mean(D) - para$fd
# mean(Y) - para$fy
# fit1 <- glm(D~X+G+Y,family = binomial(link = "logit"))
# # W0 <- matrix(diag(1/sqrt(MAFs * (1-MAFs))),ncol = nSNP0)
# # W0 <- matrix(diag(dbeta(MAF0s_true, 1, 25)), ncol = nSNP0)
# summary(fit1)
# mean(D-fit1$fitted.values) # little bias
# # round(fit1$coefficients-c(alpha0,alphaX,rep(0,ncol(G)),alphaY),3)
# fit2 <- glm(Y~X+G,family = "binomial")
# summary(fit2)
# mean(Y-fit2$fitted.values)
# round(fit2$coefficients- c(beta0,betaX,rep(0,ncol(G))),3)


GetJoint_minp <- function(D,Y,X,G,fd,fy,K=c(seq(0,0.99999,0.1),0.99)){
  n <- length(D)
  p <- ncol(X)
  Xall <- cbind(rep(1,n),X) # include intercept
  m <- ncol(G)
  MAFs <- colMeans(G)/2  # MLE of MAFs
  # W <- matrix(diag(1/sqrt(MAFs * (1-MAFs))),ncol = m)
  W <- diag(wuweights(MAFs))
  theta_hat <- GetJoint_estimation(D,Y,X,G,fd,fy)$my_est
  alphaX <- theta_hat[1:(p+1)]
  betaX <- theta_hat[(p+2):(length(theta_hat)-1)]
  alphaY <- theta_hat[length(theta_hat)]
  
  mu_D <- mu_Y <- mu_DY <- numeric(length = n)
  h11 <- h10 <- h01 <- h00 <- sum_h <- numeric(length = n)
  # h11 <- as.vector(exp(log((1-fd)/fd) + Xall%*%alphaX + alphaY + Xall%*%betaX))
  # h10 <- as.vector(exp(log((1-fd)/fd) + Xall%*%alphaX))
  
  h11 <- as.vector( exp( log( sum(D==1)*(1-fd) / (sum(D==0)*fd) ) + Xall%*%alphaX + alphaY + Xall%*%betaX) )
  h10 <- as.vector(exp( log(sum(D==1)*(1-fd)/(sum(D==0)*fd) ) + Xall%*%alphaX))
  h01 <- as.vector(exp(Xall%*%betaX))
  h00 <- rep(1,n)
  sum_h <- h11+h10+h01+h00
  
  mu_D <- as.vector((h11+h10)/sum_h)
  mu_Y <- as.vector((h11+h01)/sum_h)
  mu_DY <- as.vector(h11/sum_h)
  mu_YY <- as.vector(plogis(Xall%*%betaX))
  
  Eg <- MAFs*2
  mu_G <- matrix(rep(Eg,n),ncol = length(Eg),nrow =n, byrow = T)
  
  S_P <- c(W%*%t(G)%*%(D-mu_D),W%*%t(G)%*%(Y-mu_Y))
  S_R <- c(W%*%t(G)%*%D-W%*%t(mu_G)%*%mu_D,W%*%t(G)%*%(Y-mu_YY)-W%*%t(mu_G)%*%(mu_Y-mu_YY))
  
  
  # Phi_P
  Z14 <-  t(Xall)%*%diag(mu_D*(1-mu_D))%*%G%*%W
  Z15 <-  t(Xall)%*%diag(mu_DY-mu_D*mu_Y)%*%G%*%W
  
  Z24 <- Z15
  Z25 <- t(Xall)%*%diag(mu_Y*(1-mu_Y))%*%G%*%W
  
  Z34 <- t(mu_DY)%*%diag(1-mu_D)%*%G%*%W
  Z35 <- t(mu_DY)%*%diag(1-mu_Y)%*%G%*%W
  
  Z44 <- W%*%t(G)%*%diag(mu_D*(1-mu_D))%*%G%*%W
  # W%*%t(G)%*%diag((D-mu_D)*(D-mu_D))%*%G%*%W
  # W%*%t(mu_G)%*%diag(mu_D*(1-mu_D))%*%mu_G%*%W
  
  Z45 <-  W%*%t(G)%*%diag(mu_DY-mu_D*mu_Y)%*%G%*%W
  Z54 <- t(Z45)
  Z55 <- W%*%t(G)%*%diag(mu_Y*(1-mu_Y))%*%G%*%W
  
  ZB0 <- rbind(cbind(Z44,Z45),cbind(Z54,Z55))
  ZB1 <- rbind(cbind(Z14,Z15),cbind(Z24,Z25),cbind(Z34,Z35))
  ZB2 <- score_theta0(theta_hat,D,Y,Xall,fd)$A
  eig_ZB2 <- eigen(ZB2,symmetric = T)
  ZB2_inv <- eig_ZB2$vectors%*%diag(1/eig_ZB2$values)%*%t(eig_ZB2$vectors)
  Phi_P <- ZB0-t(ZB1)%*%ZB2_inv%*%ZB1
  
  pros <- montecarlo(score_alphaG = S_P[1:m],score_betaG = S_P[(m+1):(2*m)],
                     Phi = Phi_P,K=c(seq(0,0.99999,0.1),0.99))
  
  # Covg <- cov(G)
  # Covg <- t(G)%*%G/n - Eg%*%t(Eg)
  Covg <- var(G)
  V <- diag( 1/(MAFs*(1-MAFs)) )
  A44 <- n*V%*%Covg%*%V
  
  B2_inv <- rbind(cbind(ZB2_inv,matrix(0,nrow = (p+1)*2+1,ncol = m)),
                  cbind(matrix(0,nrow = m,ncol = (p+1)*2+1),solve(A44)))
  A15 <- t(Xall)%*%diag(mu_D*(1-mu_D))%*%mu_G%*%W
  A16 <- t(Xall)%*%diag(mu_DY-mu_D*mu_Y)%*%mu_G%*%W
  A25 <- A16
  A26 <- t(Xall)%*%diag(mu_Y*(1-mu_Y))%*%mu_G%*%W
  A35 <- t(mu_DY)%*%diag(1-mu_D)%*%mu_G%*%W
  A36 <- t(mu_DY)%*%diag(1-mu_Y)%*%mu_G%*%W
  A45 <- sum(mu_D)*V%*%Covg%*%W
  A46 <- sum(mu_Y-mu_YY)*V%*%Covg%*%W
  
  A55 <- W%*%(sum(mu_D)*Covg+sum(mu_D*(1-mu_D))*Eg%*%t(Eg))%*%W
  A56 <- W%*%(sum(mu_DY-mu_D*mu_YY)*Covg+ sum(mu_DY-mu_D*mu_Y)*Eg%*%t(Eg))%*%W
  A65 <- t(A56)
  A66 <- W%*%(sum( mu_Y*(1-mu_Y)+(mu_YY-mu_Y)^2 )*Covg+sum(mu_Y*(1-mu_Y))*Eg%*%t(Eg) )%*%W
  
  B0 <- rbind(cbind(A55,A56),cbind(A65,A66))
  B1 <- rbind(cbind(A15,A16),cbind(A25,A26),cbind(A35,A36),cbind(A45,A46))
  Phi_R <- B0-t(B1)%*%B2_inv%*%B1
  
  retr <- montecarlo(score_alphaG = S_R[1:m],score_betaG = S_R[(m+1):(2*m)],
                     Phi = Phi_R,K=c(seq(0,0.99999,0.1),0.99))
  
  #########calculate joint variance of S_P and S_R
  ##### Sigma_{XX} c(t(Xall)%*%(D-mu_D),t(Xall)%*%(Y-mu_Y),sum(D*Y)-sum(mu_DY))
  F11 <- t(Xall)%*%diag((D-mu_D)*(D-mu_D))%*%Xall-t(Xall)%*%(D-mu_D)%*%t(D-mu_D)%*%Xall/n
  F12 <- t(Xall)%*%diag((D-mu_D)*(Y-mu_Y))%*%Xall - t(Xall)%*%(D-mu_D)%*%t(Y-mu_Y)%*%Xall/n
  F13 <- t(Xall)%*%((D-mu_D)*(D*Y-mu_DY))-t(Xall)%*%(D-mu_D)%*%(sum(D*Y-mu_DY))/n
  F14 <- t(Xall)%*%diag(D-mu_D)%*%(G-mu_G)%*%V-t(Xall)%*%(D-mu_D)%*%t(rep(1,n))%*%(G-mu_G)%*%V/n
  F21 <- t(F12)
  F22 <- t(Xall)%*%diag((Y-mu_Y)*(Y-mu_Y))%*%Xall-t(Xall)%*%(Y-mu_Y)%*%t(Y-mu_Y)%*%Xall/n
  F23 <- t(Xall)%*%((Y-mu_Y)*(D*Y-mu_DY)) -t(Xall)%*%(Y-mu_Y)%*%(sum(D*Y-mu_DY))/n
  F24 <- t(Xall)%*%diag(Y-mu_Y)%*%(G-mu_G)%*%V -t(Xall)%*%(Y-mu_Y)%*%t(rep(1,n))%*%(G-mu_G)%*%V/n
  F31 <- t(F13)
  F32 <- t(F23)
  F33 <- sum(D*Y)-2*sum(mu_DY*D*Y)+sum(mu_DY^2) - sum((D*Y-mu_DY)^2)/n
  F34 <- t(D*Y-mu_DY)%*%(G-mu_G)%*%V-sum(D*Y-mu_DY)*t(rep(1,n))%*%(G-mu_G)%*%V/n
  
  M1 <- rbind(cbind(F11,F12,F13,F14),cbind(F21,F22,F23,F24),cbind(F31,F32,F33,F34))
  
  ##### Sigma_{XG}
  F15 <- t(Xall)%*%diag((D-mu_D)*D)%*%G%*%W-t(Xall)%*%diag((D-mu_D)*mu_D)%*%mu_G%*%W - t(Xall)%*%(D-mu_D)%*%t(S_R[1:m])/n
  F16 <- t(Xall)%*%diag((D-mu_D)*(Y-mu_YY))%*%G%*%W-t(Xall)%*%diag((D-mu_D)*(mu_Y-mu_YY))%*%mu_G%*%W- t(Xall)%*%(D-mu_D)%*%t(S_R[(1+m):(2*m)])/n
  F25 <- t(Xall)%*%diag((Y-mu_Y)*D)%*%G%*%W-t(Xall)%*%diag((Y-mu_Y)*mu_D)%*%mu_G%*%W- t(Xall)%*%(Y-mu_Y)%*%t(S_R[1:m])/n
  F26 <- t(Xall)%*%diag((Y-mu_Y)*(Y-mu_YY))%*%G%*%W-t(Xall)%*%diag((Y-mu_Y)*(mu_Y-mu_YY))%*%mu_G%*%W- t(Xall)%*%(Y-mu_Y)%*%t(S_R[(1+m):(2*m)])/n
  F35 <- t((D*Y-mu_DY)*D)%*%G%*%W-t((D*Y-mu_DY)*mu_D)%*%mu_G%*%W-sum(D*Y-mu_DY)*t(S_R[1:m])/n
  F36 <- t((D*Y-mu_DY)*(Y-mu_YY))%*%G%*%W-t((D*Y-mu_DY)*(mu_Y-mu_YY))%*%mu_G%*%W-sum(D*Y-mu_DY)*t(S_R[(1+m):(2*m)])/n
  
  M2 <- rbind(cbind(F15,F16),cbind(F25,F26),cbind(F35,F36))
  
  ##### Sigma_{GX}
  F41 <- W%*%t(G)%*%diag((D-mu_D)*(D-mu_D))%*%Xall-S_P[1:m]%*%t(D-mu_D)%*%Xall/n
  F42 <- W%*%t(G)%*%diag((D-mu_D)*(Y-mu_Y))%*%Xall-S_P[1:m]%*%t(Y-mu_Y)%*%Xall/n
  F43 <- W%*%t(G)%*%((D-mu_D)*(D*Y-mu_DY))-sum(D*Y-mu_DY)*S_P[1:m]/n
  F44 <- W%*%t(G)%*%diag(D-mu_D)%*%(G-mu_G)%*%V-S_P[1:m]%*%t(rep(1,n))%*%(G-mu_G)%*%V/n
  F51 <- W%*%t(G)%*%diag((D-mu_D)*(Y-mu_Y))%*%Xall-S_P[(1+m):(2*m)]%*%t(D-mu_D)%*%Xall/n
  F52 <- W%*%t(G)%*%diag((Y-mu_Y)*(Y-mu_Y))%*%Xall-S_P[(1+m):(2*m)]%*%t(Y-mu_Y)%*%Xall/n
  F53 <- W%*%t(G)%*%((Y-mu_Y)*(D*Y-mu_DY))-sum(D*Y-mu_DY)*S_P[(1+m):(2*m)]/n
  F54 <- W%*%t(G)%*%diag(Y-mu_Y)%*%(G-mu_G)%*%V-S_P[(1+m):(2*m)]%*%t(rep(1,n))%*%(G-mu_G)%*%V/n
  
  M3 <- rbind(cbind(F41,F42,F43,F44),cbind(F51,F52,F53,F54))
  
  ##### Sigma_{GG}
  
  F45 <- W%*%t(G)%*%diag((D-mu_D)*D)%*%G%*%W-W%*%t(G)%*%diag((D-mu_D)*mu_D)%*%mu_G%*%W-S_P[1:m]%*%t(S_R[1:m])/n
  F46 <- W%*%t(G)%*%diag((D-mu_D)*(Y-mu_YY))%*%G%*%W-W%*%t(G)%*%diag((D-mu_D)*(mu_Y-mu_YY))%*%mu_G%*%W-S_P[1:m]%*%t(S_R[(1+m):(2*m)])/n
  F55 <- W%*%t(G)%*%diag((Y-mu_Y)*D)%*%G%*%W-W%*%t(G)%*%diag((Y-mu_Y)*mu_D)%*%mu_G%*%W-S_P[(1+m):(2*m)]%*%t(S_R[1:m])/n
  F56 <- W%*%t(G)%*%diag((Y-mu_Y)*(Y-mu_YY))%*%G%*%W-W%*%t(G)%*%diag((Y-mu_Y)*(mu_Y-mu_YY))%*%mu_G%*%W-S_P[(1+m):(2*m)]%*%t(S_R[(1+m):(2*m)])/n
  
  M4 <- rbind( cbind(F45,F46),cbind(F55,F56))
  star_Sigma_PP <- Phi_P
  star_Sigma_PR <- M4 - t(ZB1)%*%ZB2_inv%*%M2 - M3%*%B2_inv%*%B1 + t(ZB1)%*%ZB2_inv%*%M1%*%B2_inv%*%B1
  # star_Sigma_PR <- Sigma_PR-cbind(t(ZB1),Sigma_PM)%*%B2_inv%*%B1
  star_Sigma_RR <- Phi_R
  
  joint_sigma <- rbind(cbind(star_Sigma_PP,star_Sigma_PR),cbind(t(star_Sigma_PR),star_Sigma_RR))
  
  ############calculate S_EB1 and Phi_EB1
  psi <- S_P-S_R
  K1 <- (psi%*%t(psi))%*%solve(psi%*%t(psi)+Phi_P)
  S_EB1 <- K1%*%S_P + (diag(1,nrow = 2*m)-K1)%*%S_R
  deno <- 1 + as.vector(t(psi)%*%solve(Phi_P)%*%psi)
  dd1 <- diag(rep(1/deno,2*m))-2*psi%*%t(psi)%*%solve(Phi_P)/(deno^2)
  # dd <- diag( (diag(Phi_P)*(diag(Phi_P)-psi^2)) / ((diag(Phi_P)+psi^2)^2) )
  cc1 <- cbind( diag(1,nrow = 2*m)-dd1,dd1) 
  Phi_EB1 <- cc1%*%joint_sigma%*%t(cc1)
  EB1 <- montecarlo(score_alphaG = S_EB1[1:m],score_betaG = S_EB1[(m+1):(2*m)],
                    Phi = Phi_EB1,K=c(seq(0,0.99999,0.1),0.99))

  ############## Two step
  p_ind <- GetInd_pval(D,G,X,1000)
  if(p_ind < 0.05/2){
    Ts1 <- pros
  }else{
    Ts1 <- retr
  }
  
  if(p_ind < 0.1){
    Ts2 <- pros
  }else{
    Ts2 <- retr
  }

  return(list(pros=pros,retr=retr,EB=EB1,Ts1 = Ts1,Ts2 = Ts2))
}






montecarlo <- function(score_alphaG,score_betaG,Phi,K=c(seq(0,0.99999,0.1),0.99)){
  m <- length(score_alphaG)
  Phi11 <- Phi[(1:m),(1:m)]
  Phi12 <- Phi[(1:m),(m+1):(2*m)]
  Phi21 <- Phi[(m+1):(2*m),(1:m)]
  Phi22 <- Phi[(m+1):(2*m),(m+1):(2*m)]
  Qk <- p <- numeric(length = length(K))
  lambda <- matrix(NA,nrow = length(K),ncol = 2*m)
  for(i in 1:length(K)){ # compute Q_k for each k 
    k <- K[i]
    Qk[i] <- k * t(score_alphaG)%*%score_alphaG +(1-k)*t(score_betaG)%*%score_betaG
    
    Sigma_k <- diag(c(rep(sqrt(k),m),rep(sqrt(1-k),m)))%*%Phi%*%diag(c(rep(sqrt(k),m),rep(sqrt(1-k),m)))
    # Sigma_k <- rbind(cbind(k*Phi11,sqrt(k*(1-k))*Phi12),cbind(sqrt(k*(1-k))*Phi21,(1-k)*Phi22))
    lambda[i,] <- eigen(Sigma_k,symmetric = TRUE)$values
    if(any(lambda[i,] > 0)){
      p[i] <- pchisqsum2(Qk[i],lambda[i,],method = "saddlepoint")$p
    }else{
      p[i] <- 1
    }
  }
  # Minimum p-value test
  Tobs <- min(p)
  kopt <- K[which.min(p)]
  q <-  upper <- numeric(length = length(K))
  for(i in 1:length(K)){ # calculate (1-Tobs) quantiles of Q_k
    sat <- satterthwaite(lambda[i,],rep(1,length(lambda[i,])))
    upper[i] <- stats::qchisq(Tobs/20,df=sat$df,lower.tail = FALSE)*sat$scale*1000
    tmpT <- try(stats::uniroot(function(x){pchisqsum2(x,lambda=lambda[i,],method="liu",acc=1e-7)$p- Tobs }, interval=c(1e-10,upper[i]))$root, silent = TRUE)
    if(class(tmpT) == "try-error"){
      q[i] <- Qk[i]
      errflag <- 2
    }else{
      q[i] <- tmpT
    }
    # q[i] <- Getqk(lambda[i,],Tobs,upper)
  }
  Phi11_inv <- solve(Phi11)
  eig <- eigen(Phi22-Phi21%*%Phi11_inv%*%Phi12,symmetric = T)
  lambda <- eig$values
  D <- diag(lambda)
  diag(D)[zapsmall(diag(D)) > 0] <- 1/sqrt(diag(D)[zapsmall(diag(D)) > 0])
  diag(D)[diag(D) <= 0] <- 0
  Uvec <- eig$vectors
  Udel <- D%*%t(Uvec)%*%Phi21%*%Phi11_inv
  
  Fcond <- function(x){
    qmeta <- min((q-K*sum(x^2))/(1-K))
    delta <- c(Udel%*%x)^2
    return(pchisqsum2(qmeta,lambda,delta,method = "saddlepoint")$p)
  }
  
  B <- 1e4
  # set.seed(1234)
  xx <- MASS::mvrnorm(B,mu = rep(0,m),Phi11)
  pmin <- apply(xx, 1, Fcond)
  
  return(list(pactual = mean(pmin, na.rm = T), min_p = Tobs,kopt = kopt))
}

GetJoint_estimation <- function(D,Y,X,G,fd,fy){
  n <- length(D)
  p <- ncol(X)
  Xall <- cbind(rep(1,n),X) # include intercept
  theta0 <- rep(0,(p+1)*2+1)  # Parameters to be estimated
  
  # initial value
  fit1 <- glm(D~X+Y,family = "binomial")
  fit2 <- glm(Y~X,family = "binomial")
  alphaX0 <- c(fit1$coefficients[1]-log(sum(D==1)*(1-fd)/(sum(D==0)*fd)),fit1$coefficients[-c(1,p+2)])
  betaX0 <- c(fit2$coefficients[1]-log(sum(Y==1)*(1-fy)/(sum(Y==0)*fy)),fit2$coefficients[-1])
  alphaY0 <- fit1$coefficients[p+2]
  
  theta0 <- as.vector(c(alphaX0,betaX0,alphaY0))
  # Sigma0 <- as.vector(c(summary(fit1)$coefficients[,"Std. Error"][-(p+2)],
  #                       summary(fit2)$coefficients[,"Std. Error"],
  #                       summary(fit1)$coefficients[,"Std. Error"][(p+2)])) #sd
  v1 <- diag(vcov(fit1))
  v2 <- diag(vcov(fit2))
  Sd0 <- sqrt(c(v1[-(p+2)],v2,v1[(p+2)])) #cov
  Newtons <- function(theta0,eps = 1e-11,it_max=100){
    index <- 0
    k <- 1
    
    while(k<=it_max){
      
      score0 <- score_theta0(theta0,D,Y,Xall,fd)
      theta1 <- theta0 + solve(score0$A)%*%score0$score
      norm <- norm(theta1-theta0,type = "2")
      if(norm < eps){
        index <- 1 ;break 
      }
      k <- k+1
      theta0 <- theta1
      
    }
    list(root=theta1, it = k ,index=index)
  }
  # Newtons(theta0)
  theta_hat <- Newtons(theta0)$root
  Sd_hat <- sqrt(diag(solve(score_theta0(theta_hat,D,Y,Xall,fd)$A)))
  
  return(list(glm_est = as.vector(c(alphaX0,betaX0,alphaY0)),glm_sd =Sd0, 
              my_est = theta_hat,my_sd=Sd_hat))
}


satterthwaite <- function(a, df){
  if(any(df > 1)){
    a <- rep(a,df)
  }
  tr <- mean(a)
  tr2 <- mean(a^2)/(tr^2)
  list(scale = tr * tr2,df = length(a)/tr2)
}

wuweights <- function(maf){
  ifelse(maf>0, dbeta(maf, 1, 25), 0)
} 

# calculate (1-Tobs) quantiles of Q_k
# Getqk <- function(lambda,Tobs,upper){
#   fn <- function(x){
#     if(any(lambda > 0)){
#       p <- pchisqsum2(x,lambda,method = "liu",acc=1e-5)$p
#     }else{
#       p <- 1
#     }
#     return(p-Tobs)
#   }
#   qk <- uniroot(fn,c(1e-10,upper=upper))$root
#   return(qk)
# }

  
# calculate score function and Fisher Information under H0
score_theta0 <- function(theta,D,Y,Xall,fd){ #p+1,p+1,1
  n <- length(D)
  p <- ncol(Xall)-1
  mu_D <- mu_Y <- mu_DY <- numeric(length = n)
  h11 <- h10 <- h01 <- h00 <- sum_h <- numeric(length = n)
  alphaX <- theta[1:(p+1)]
  betaX <- theta[(p+2):(length(theta)-1)]
  alphaY <- theta[length(theta)]
  # print(dim(Xall))
  # print(alphaX)
  # print(Xall%*%alphaX)
  
  h11 <- as.vector( exp( log( sum(D==1)*(1-fd) / (sum(D==0)*fd) ) + Xall%*%alphaX + alphaY + Xall%*%betaX) )
  h10 <- as.vector(exp( log(sum(D==1)*(1-fd)/(sum(D==0)*fd) ) + Xall%*%alphaX))
  h01 <- as.vector(exp(Xall%*%betaX))
  h00 <- rep(1,n)
  sum_h <- h11+h10+h01+h00

  mu_D <- as.vector((h11+h10)/sum_h)
  mu_Y <- as.vector((h11+h01)/sum_h)
  mu_DY <- as.vector(h11/sum_h)
  
  # mean(D-mu_D)
  # mean(Y-mu_Y)
  # mean(D*Y-mu_DY)

  score <- c(t(Xall)%*%(D-mu_D),t(Xall)%*%(Y-mu_Y),sum(D*Y)-sum(mu_DY))
  
  # Fisher Information
  A11 <- t(Xall)%*%diag(mu_D*(1-mu_D))%*%Xall
  A12 <- t(Xall)%*%diag(mu_DY-mu_D*mu_Y)%*%Xall
  A13 <- t(Xall)%*%diag(1-mu_D)%*%mu_DY

  A21 <- t(A12)
  A22 <- t(Xall)%*%diag(mu_Y*(1-mu_Y))%*%Xall
  A23 <- t(Xall)%*%diag(1-mu_Y)%*%mu_DY

  A31 <- t(A13)
  A32 <- t(A23)
  A33 <- sum(mu_DY*(1-mu_DY))

  A <- rbind(cbind(A11,A12,A13),cbind(A21,A22,A23),cbind(A31,A32,A33)) # Sigma
  
  return(list(score = score, A = A)) #, mu_D = mu_D, mu_Y = mu_Y ,mu_DY = mu_DY
}  
  


