library(data.table)
source("joint_test.R")
Pheno_covariate_data <- fread(file = "PretermData/PD_QC_MCpairs_phenotype.txt", header = TRUE, sep = " ")  ## 1626 samples have maternal information
BMINA_ID <- which(is.na(Pheno_covariate_data$maternal_bmi))
pheno_1 <- Pheno_covariate_data[-BMINA_ID, ] # 1520*52
BirthwtNA_ID <- which(is.na(pheno_1$birthweight))
pheno_deleteNA <- pheno_1[-BirthwtNA_ID, ] # 1508*52

# primary phenotype
D <- pheno_deleteNA$premie; table(D) #851control,657case
# secondary phenotype
Y <- ifelse(pheno_deleteNA$birthweight<3280,1,0);table(Y) # 309case,1199control
#Y <- ifelse(pheno_deleteNA$birthweight<2500,1,0);table(Y)# 2500
# environmental factors
BMI <- scale(pheno_deleteNA$maternal_bmi); # 这里的环境因素或者次级表型能不能选吸烟?,没有na
# covariates
m_height <- scale(pheno_deleteNA$maternal_height);  # mean = 169
m_age <- scale(pheno_deleteNA$maternal_age)  # mean = 30
# m_height <- pheno_deleteNA$maternal_height;  # mean = 169
# m_age <- pheno_deleteNA$maternal_age  # mean = 30
load(file="~/PretermData/MGdata.RData")

birthwtG <- MGdata[-BirthwtNA_ID,]

X <- cbind(m_age,m_height,BMI)
dist <- 59275
TEKT3_SNPs <- intersect(which(position >= (15207129-dist) & position <= (15244914-dist)) , which(chro==17))
TEKT3 <- all_SNP[TEKT3_SNPs]
G <- as.matrix(birthwtG)[, TEKT3_SNPs]

dist <- 2182477 # grch37-本数据=2182477
IGF1R_SNPs <- intersect(which(position >= (99191768-dist) & position <= (99507759-dist)) , which(chro==15))
IGF1R <- all_SNP[IGF1R_SNPs]
G <- as.matrix(birthwtG)[, IGF1R_SNPs]

dist <- 127413
WNT4_SNPs <- intersect(which(position >= (22443806-dist) & position <= (22469590-dist)) , which(chro==1))
WNT4 <- all_SNP[WNT4_SNPs]
G <- as.matrix(birthwtG)[, WNT4_SNPs]

pvalue <- numeric(4)
p_ind <- GetInd_pval(D,G,X,1000)
pvalue[1] <- p_ind
# GetJoint_estimation(D,Y,X,G,0.05,0.1)

my <- GetJoint_minp(D,Y,X,G,fd = 0.05,fy = 0.2,K=c(seq(0,0.99999,0.05),0.99))
pvalue[2] <- my$pros$pactual
pvalue[3] <- my$retr$pactual
pvalue[4] <- my$Ts$pactual