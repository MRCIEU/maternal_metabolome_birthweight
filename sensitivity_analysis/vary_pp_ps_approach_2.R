ptm=proc.time()

library(readr)
library(ggplot2)
library(data.table)
library(TwoSampleMR)
library(stringr)

source("VZ_summary_mvMR_SSS_function.R")
source("VZ_summary_mvMR_BF_function.R")

ao=available_outcomes()
mvmr_harm=read_csv("mvmr_harm.csv")
nmr_metabolites_UKBB=read_csv("nmr_metabolites.csv")

#reshape dataframe
setDT(mvmr_harm)
test=reshape(mvmr_harm, timevar="id.exposure", idvar=c("SNP"), direction="wide")
test=as.matrix(test)
col=as.numeric(grep("exposure|SNP", colnames(test)))
test=test[,col]
length(unique(test[,1]))
snps=unique(mvmr_harm$SNP)
effect_alleles=as.data.frame(matrix(1:length(snps),nrow=length(snps),ncol=3))
colnames(effect_alleles)=c("SNP", "effect_allele.exposure", "other_allele.exposure")
for (i in 1:length(snps))
{
  effect_alleles[i,]=mvmr_instruments_ukbb[i, c("SNP", "effect_allele.exposure", "other_allele.exposure")]
}
merged=merge(effect_alleles, test, by="SNP")
outcome_dat=read_outcome_data(snps=merged$SNP, filename="Maternal_Effect_European_meta_NG2019.txt", beta_col = "beta", effect_allele_col = "ea", other_allele_col = "nea" ,
                               eaf_col = "eaf", pval_col = "p", phenotype_col = "Birthweight")
rs=merged$SNP
bw_beta=outcome_dat$beta.outcome
bw_se=outcome_dat$se.outcome
find_beta_cols=grep("beta.exposure",names(merged))
betaX_2=merged[,find_beta_cols]
rf_2=colnames(betaX_2)
rs_2=merged[,1]
betaX_2.2=sapply(betaX_2, as.numeric )        
betaX_ivw_2=betaX_2.2/bw_se
bw_beta_ivw_2=bw_beta/bw_se
bw_nmr_input_2=new(Class = "mvMRInput", betaX = as.matrix(betaX_ivw_2), betaY=as.matrix(bw_beta_ivw_2), snps=rs_2, exposure=rf_2, outcome = "Birthweight")

ao=available_outcomes()
ao=ao[grep("met-d", ao$id),]

prior_prob=c(0.01, 0.05, 0.10, 0.20, 0.30)
mr.bw_BMA.out.step1=vector("list",length(prior_prob))
best.model.out.step1=vector("list",length(prior_prob))
for (i in 1:length(prior_prob))
{
  BMA_output=summarymvMR_SSS(bw_nmr_input_2, kmin=1, kmax=12, prior_prob=prior_prob[i], max_iter=100000)
  mr.bw_BMA.out.step1[[i]]=sss.report.mr.bma(BMA_output, top=10)
  best.model.out.step1[[i]]=sss.report.best.model(BMA_output, prior_sigma=0.5, top=10)
  print(i)
}
for (i in 1:length(best.model.out.step1))
{
  rename_traits=as.data.frame(best.model.out.step1[[i]])
  rename_traits$`rf combination`=gsub("beta.exposure_","",rename_traits$`rf combination`)
  rename_traits_mr=as.data.frame(mr.bw_BMA.out.step1[[i]])
  rename_traits_mr$`rf combination`=gsub("beta.exposure_","",rename_traits$`rf combination`)
  for (j in 1:nrow(ao))
  {
    rename_traits$`rf combination`=gsub(ao$id[j], ao$trait[j], rename_traits$`rf combination`)
    rename_traits_mr$`rf combination`=gsub(ao$id[j],ao$trait[j], rename_traits_mr$`rf combination`)
    rename_traits$`rf combination`=gsub("beta.exposure.","", rename_traits$`rf combination`)
    rename_traits_mr$`rf combination`=gsub("beta.exposure.","", rename_traits_mr$`rf combination`)
  }
  assign(paste0("best_model_out_step1",i),rename_traits)
  assign(paste0("mr_bw_BMA.out_step1",i),rename_traits)
}
best.model.out.step1.updated=do.call("list",mget(grep("best_model_out_step1",names(.GlobalEnv),value=TRUE)))
mr.bw_BMA.out.step1.updated=do.call("list",mget(grep("mr_bw_BMA.out_step1",names(.GlobalEnv),value=TRUE)))

sigma=c(0.1, 0.3, 0.5, 0.7)
mr.bw_BMA.out.step3=vector("list",length(sigma))
best.model.out.step3=vector("list",length(sigma))
for (i in 1:length(sigma))
{
  BMA_output=summarymvMR_SSS(bw_nmr_input_2, kmin=1, kmax=12, prior_prob=0.1, max_iter=100000)
  mr.bw_BMA.out.step3[[i]]=sss.report.mr.bma(BMA_output, top=10)
  best.model.out.step3[[i]]=sss.report.best.model(BMA_output, prior_sigma=sigma[i], top=10)
  print(i)
}
for (i in 1:length(best.model.out.step3))
{
  rename_traits=as.data.frame(best.model.out.step3[[i]])
  rename_traits$`rf combination`=gsub("beta.exposure_","",rename_traits$`rf combination`)
  rename_traits_mr=as.data.frame(mr.bw_BMA.out.step3[[i]])
  rename_traits_mr$`rf combination`=gsub("beta.exposure_","",rename_traits$`rf combination`)
  for (j in 1:nrow(ao))
  {
    rename_traits$`rf combination`=gsub(ao$id[j], ao$trait[j], rename_traits$`rf combination`)
    rename_traits_mr$`rf combination`=gsub(ao$id[j], ao$trait[j], rename_traits_mr$`rf combination`)
    rename_traits$`rf combination`=gsub("beta.exposure.","", rename_traits$`rf combination`)
    rename_traits_mr$`rf combination`=gsub("beta.exposure.","", rename_traits_mr$`rf combination`)
  }
  assign(paste0("best_model_out_step3",i),rename_traits)
  assign(paste0("mr_bw_BMA.out_step3",i),rename_traits)
}
best.model.out.step3.updated=do.call("list",mget(grep("best_model_out_step3",names(.GlobalEnv),value=TRUE)))
mr.bw_BMA.out.step3.updated=do.call("list",mget(grep("mr_bw_BMA.out_step3",names(.GlobalEnv),value=TRUE)))

superset=do.call(cbind, Map(cbind, best.model.out.step1.updated))
df=data.frame(matrix(unlist(superset),nrow=10))
RF=seq(1,15,by=3)
PP=seq(2,16,by=3)
CE=seq(3,17,by=3)
colnames(df)[RF]=paste0("RF_model",1:5,"")
colnames(df)[PP]=paste0("Posteriorprob_model",1:5)
colnames(df)[CE]=paste0("Causalestimate_model",1:5)

superset=do.call(cbind, Map(cbind, mr.bw_BMA.out.step1.updated))
df=data.frame(matrix(unlist(superset),nrow=10))
RF=seq(1,15,by=3)
MI=seq(2,16,by=3)
AE=seq(3,17,by=3)
colnames(df)[RF]=paste0("RF_model",1:5,"")
colnames(df)[MI]=paste0("MarginalInclusion_model",1:5)
colnames(df)[AE]=paste0("Averageeffect_model",1:5)

superset=do.call(cbind, Map(cbind, best.model.out.step3.updated))
df=data.frame(matrix(unlist(superset),nrow=10))
RF=seq(1,12,by=3)
PP=seq(2,13,by=3)
CE=seq(3,14,by=3)
colnames(df)[RF]=paste0("RF_model",1:4,"")
colnames(df)[PP]=paste0("Posteriorprob_model",1:4)
colnames(df)[CE]=paste0("Causalestimate_model",1:4)

superset=do.call(cbind, Map(cbind, mr.bw_BMA.out.step3.updated))
df=data.frame(matrix(unlist(superset),nrow=10))
RF=seq(1,12,by=3)
MI=seq(2,13,by=3)
AE=seq(3,14,by=3)
colnames(df)[RF]=paste0("RF_model",1:4,"")
colnames(df)[MI]=paste0("MarginalInclusion_model",1:4)
colnames(df)[AE]=paste0("Averageeffect_model",1:4)

# model diagnostics: main analyses - prior sigma = 0.5, prior prob = 0.1, best.model.out
diag_ppthresh=0.0004
best.model.out=best.model.out.step3.updated[[3]]
nr_diag=length(which(best.model.out[,2]>=diag_ppthresh))
model_index=names(which(best.model.out[,2]>=diag_ppthresh))
title=rep("1", nr_diag)
predicted_bw=matrix(ncol=nr_diag, nrow=length(bw_beta_ivw_2))
cD=matrix(ncol=nr_diag, nrow=length(bw_beta_ivw_2))
cD_thresh=0
Q=matrix(ncol=nr_diag, nrow=length(bw_beta_ivw_2))

for(i in 1:nr_diag){
  print(as.numeric(unlist(strsplit(model_index[i], ","))))
  if(length(as.numeric(unlist(strsplit(model_index[i], ","))))>1){
    betaX_model=cbind(betaX_ivw_2[,as.numeric(unlist(strsplit(model_index[i], ",")))])
  }
  else{
    betaX_model=as.matrix(betaX_ivw_2[,as.numeric(unlist(strsplit(model_index[i], ",")))])
  }	
  title[i]=paste(rf_2[as.numeric(unlist(strsplit(model_index[i], ",")))],collapse=' + ')
  sigma_vec=rep(0.5, ncol(betaX_model))
  cD[,i]=cooksD(bw_beta_ivw_2, betaX_model,sigma_vec)$cooksD
  cD_thresh[i]=cooksD(bw_beta_ivw_2, betaX_model, sigma_vec)$cooksD_thresh
  H_fm=betaX_model%*% solve(t(betaX_model) %*% betaX_model + sigma_vec^{-2} ) %*% t(betaX_model)
  predicted_bw[,i]=H_fm %*% bw_beta_ivw_2
  Q[,i] = (bw_beta_ivw_2-predicted_bw[,i])^2
}
maxCD=apply(cD, MARGIN=1, FUN=max)
sort.ix=sort.int(maxCD, decreasing=TRUE, index.return=TRUE)
cooksD_tab=cbind(rs_2,round(cD,digits=3))
cooksD_top_30=as.data.frame(cooksD_tab[sort.ix$ix,][1:30,])
colnames(cooksD_top_30)=c("rs","cooksd")

cD_thresh
for(i in 1:nr_diag){
  print(rs_2[which(cD[,i] > cD_thresh[i])])
}

maxQ=apply(Q, MARGIN=1, FUN=max)
sort.ix = sort.int(maxQ, decreasing=TRUE, index.return=TRUE)
Q_tab=cbind(rs_2,round(Q,digits=3))
Q_stat_top_30=as.data.frame(Q_tab[sort.ix$ix,][1:30,])
colnames(Q_stat_top_30)=c("rs","qstat")

# Creating plots 
ao=available_outcomes()
trait_linkage=ao[grep("met-d", ao$id), c("id","trait")]

for (j in 1:nrow(trait_linkage))
{
  title=gsub(trait_linkage$id[j], trait_linkage$trait[j], title)
}
title=gsub("beta.exposure.","",title)

plot_list=list()
pdf("cooks_d_plots")
for(i in 1:nr_diag){
  df = data.frame(x=predicted_bw[,i], y =bw_beta_ivw_2, cD = cD[,i], rs=rs_2[which(rs_2%in%outcome_dat$SNP)])
  p=ggplot(df, aes(x, y)) + geom_point(aes(colour = cD), size =4) +
    scale_colour_gradientn(colours = c("white", "orange", "red", "darkred"))  +
    labs(x = "predicted beta bwt", y="observed beta bwt", colour="Cooks D") + geom_hline(yintercept = 0, linetype="dotted") +
    geom_vline(xintercept = 0, linetype="dotted")  +
    geom_text(aes(label=ifelse(cD>cD_thresh[i],as.character(rs[which(rs%in%outcome_dat$SNP)]),'')),hjust=0.5, vjust=-1, size=5)+
    theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13), axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.text=element_text(size=16),legend.title=element_text(size=18))+
    ggtitle(title[i])
  plot_list[[i]]=p
  print(plot_list[[i]])
}
dev.off

plot_list=list()
pdf("q_plots")
for(i in 1:nr_diag){
  
  df = data.frame(x=predicted_bw[,i], y =bw_beta_ivw_2, Q = Q[,i], rs=rs_2[which(rs_2%in%outcome_dat$SNP)])
  p = ggplot(df, aes(x , y)) +  geom_point(aes(colour = Q), size =4) + scale_colour_gradientn(colours =
                                                                                                c("#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"))  +
    labs(x = "predicted beta bwt", y="observed beta bwt", colour="Q") +
    geom_hline(yintercept = 0, linetype="dotted") + geom_vline(xintercept = 0, linetype="dotted") +
    geom_text(aes(label=ifelse(Q>10,as.character(rs[which(rs%in%outcome_dat$SNP)]),'')),hjust=0.5, vjust=-1, size=5) +
    theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13), axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.text=element_text(size=16),legend.title=element_text(size=18)) +
    ggtitle(title[i])
  plot_list[[i]]=p
  print(plot_list[[i]])
}
dev.off()

proc.time()-ptm
