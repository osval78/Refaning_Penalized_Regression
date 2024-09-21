rm(list=ls(all=TRUE))
library(dplyr)
library(writexl)
source('glmnet_1Pred_1Pen-v1.R')

Dir = '.'
Folds_df = read.csv(file=paste0(Dir,'/10Folds_df.csv'))
DataSets = c('Disease_AL','EYT_1_AL','EYT_2_AL','EYT_3_AL','Groundnut_AL',
             'Indica_AL','Japonica_AL','Maize_AL','Wheat_1_AL','Wheat_2_AL',
             'Wheat_3_AL','Wheat_4_AL','Wheat_5_AL','Wheat_6_AL')  
head(Folds_df)
Traits_df = read.csv(file=paste0(Dir,'/Traits_df.csv'))

head(Folds_df)
#DataSets =  DataSets[1:9]
Tb =  data.frame(); Preds = Tb
d = 1#Choose the dataset
load(paste0(DataSets[d],'.RData'),verbose=TRUE)
head(Pheno); Pheno$GID =  Pheno$Line; #head(Pheno)
ZL =  model.matrix(~0+GID,data=Pheno); Indicador = mean(colnames(Geno)==sub('GID','',colnames(ZL)))
XR = XR2_f(Geno)#KL = ZL%*%as.matrix(Geno)%*%t(ZL)
Fold_df = droplevels(Folds_df[Folds_df$DataSet==DataSets[d],])
head(Fold_df)
Traits_d= Traits_df$Traits[Traits_df$DataSet==DataSets[d]]
for(t in 1:length(Traits_d)){
  Pheno$Response = Pheno[,Traits_d[t]]
  for(k in 1:length(unique(Fold_df$Fold))){
    set.seed(10)  
    Pos_tst =  Fold_df$Pos_tst[Fold_df$Fold==k]; y =  Pheno$Response; y_tst = y[Pos_tst]
    Pheno$Ind_tr = TRUE; Pheno$Ind_tr[Pos_tst] = FALSE
    y_tr =  y[-Pos_tst]
    #Method: glmnet-M
    Time2 = proc.time()
    FM_ls = glmnet1_f(XR[-Pos_tst,],y[-Pos_tst])#cv.glmnet(XR[Pos_tr,],y_tr,alpha = 0,standardize=F)#,foldid =IFold)#,lambda = lambda1)#,parallel = TRUE)
    yp2_tst = predict(FM_ls$FM, newx= XR[Pos_tst,],s='lambda.min')[,1]
    Time2 = proc.time()-Time2
    Preds2_k = data.frame(Model='glmnet-M',Fold = k,y = y_tst,yp = yp2_tst,
                          Time = Time2[[3]]) 
    #Method: glmnet
    Time0 = proc.time()
    FM = cv.glmnet(XR[-Pos_tst,],y[-Pos_tst],alpha = 0,standardize=F)#,foldid = FM_ls$IFold)#,lambda = lambda1)#,parallel = TRUE)
    yp0_tst = predict(FM, newx= XR[Pos_tst,],s='lambda.min')[,1]
    #plot(y_tst,yp0_tst); abline(a=0,b=1)
    Time0 = proc.time()-Time0
    Preds0_k = data.frame(Model='glmnet',Fold = k,y = y_tst,yp = yp0_tst,
                          Time = Time0[[3]]) 
    Preds_k =  rbind(Preds0_k,Preds2_k)
    plot(yp0_tst,yp2_tst); abline(a=0,b=1)
    Preds =  rbind(Preds,data.frame(Data = DataSets[d],Trait = Traits_d[t],Preds_k))
    head(Preds)
    cat('k=',k,'\n')
  }
}

Tb =  Preds%>%group_by(Model,Trait,Fold)%>%summarise(MSE =  mean((y-yp)^2),Cor =  cor(y,yp),
                                                NRMSE = sqrt(mean((y-yp)^2))/mean(y),
                                                Time =  mean(Time))
Tb_Smm = data.frame(Tb%>%group_by(Model,Trait)%>%summarise(MSE_Mean = mean(MSE), MSE_SD =  sd(MSE),
                                              Cor_Mean =  mean(Cor), Cor_SD =  sd(Cor),
                                              NRMSE_Mean = mean(NRMSE), NRMSE_SD =  sd(NRMSE),
                                              Time_Mean =  mean(Time), Time_SD = sd(Time) ))
Tb_Smm
