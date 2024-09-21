library(splitTools)
library(glmnet)
Vars_f<-function(s2y,X,nGrid=100,  a=0.001, b=0.999)
{
  Mean_trtXX = mean(apply(X,1,function(x)sum(x^2)))# mean(diag(t(X)%*%X2)#
  R2v =  exp(seq(log(a),log(b),length=nGrid))#(b-a)/nGrid))
  sigma2Bv = (R2v)*s2y/Mean_trtXX; sigma2v = (1-R2v)*s2y
  #plot(delta2v,sigma2Bv)
  list(lambda=sigma2v/sigma2Bv,sigma2B=sigma2Bv)
}
SKFCV_f<-function(y,K=5,Reps=1,n_bins=20,seed=1)
{
  F_tst_ls = create_folds(y,k = K, type ='stratified',invert=TRUE,#c("stratified", "basic", "grouped", "blocked"),
                          shuffle = FALSE,seed = seed, n_bins = n_bins,m_rep = Reps)#, m_rep = 1,use_names = TRUE, invert = FALSE,
  #table(findInterval(y,[IF_tst_ls[[1]]],quantile(y_tr,seq(0.05,0.95,0.05))))
  Fold =  1:length(y); for(k in 1:K) Fold[F_tst_ls[[k]]] = k
  Fold
}
library(MASS)
glmnet1_f <- function(X_tr,y_tr,Params_ls=NULL,Plot=FALSE,trace.it=0)#Inputs_ls =  list(Input1, Input2, Input3)
{ #Pars =  list(nGrid1=1)
  Params = Params_ls; Params_ls = NULL; NPars =  names(Params); 
  Params0 = list('nGrid1'=100,'a1'=0.00001,'b1'=0.9999,'IK'= 10,'Reps'=1,'n_bins'=20)
  PosD = which(!(names(Params0)%in%NPars)); Params[names(Params0)[PosD]] = Params0[PosD]
  lambv = Vars_f(var(y_tr),X_tr,nGrid = Params$nGrid1,  a = Params$a1, 
                 b = Params$b1)$lambda
  IFold = SKFCV_f(y_tr,K = Params$IK,Reps = Params$Reps,n_bins = Params$n_bins)
  n_tr =  nrow(X_tr)
  FM = cv.glmnet(X_tr,y_tr,alpha = 0,trace.it = trace.it,standardize=FALSE,
                 foldid = IFold,lambda = lambv/(n_tr-sum(IFold==100)))#,parallel = TRUE)
    lambda_min = FM$lambda.min
  if(Plot==TRUE) plot(FM)
  list(FM=FM,lambda1_O = lambda_min,IFold =  IFold,lambda = lambv/n_tr)
}

XR2_f<-function(G,Pheno = NULL)
{
  svd_G =  svd(G);R = (svd_G$u%*%diag(sqrt(svd_G$d))); PosR = which(svd_G$d>1e-8)
  X_R = (svd_G$u[,PosR]%*%diag(svd_G$d[PosR]))
  if(!is.null(Pheno)){
    ZL =  model.matrix(~0+GID,data=Pheno); Pos =  match(colnames(ZL),row.names(G))
    X_R =  ZL%*%X_R[Pos,]
  }
  X_R
}