library(aSPU)
GetInd_pval <- function(D, G, X,  n.perm = n.perm){
  
  G <- as.matrix(G)
  X <- as.matrix(X)
  m <- ncol(G)
  
  MAFs <- colMeans(G, na.rm = T)/2
  E.inter <- as.vector(scale(X[, ncol(X)], center = T, scale = T))
  id <- which(D==0)
  cov <- X[id, -(ncol(X))]
  Y <- E.inter[id]
  X <- G[id,]
  ### Test the independence between E and causal SNPs 
  rest <- aSPU::aSPU(Y=Y, X=X, cov = cov, resample = "perm", model = "gaussian", n.perm = 1000)
  p0 <- rest$pvs["aSPU"]
  #pval <- apply(G[id, ], 2, function(x) kruskal.test(E.inter[id]~x)$p.value) 
  #p0 <- pcauchy(as.numeric(mean(tan((0.5-pval)*pi))), scale = 1, lower.tail = FALSE)
  return(p0)
}