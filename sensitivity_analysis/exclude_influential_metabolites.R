# excluding glycolysis related metabolites
library(readr); library(ggplot2); library(data.table);library(TwoSampleMR);library(stringr);library(readxl)
library(tidyverse); library(ggforestplot)
source("VZ_summary_mvMR_SSS_function.R")
source("VZ_summary_mvMR_BF_function.R")

## UKBB
setwd()
nmr_metabolites=read_excel("nmr_metabolites.xlsx")
colnames(nmr_metabolites)[1]="GWAS_id"
nmr_metabolites$GWAS_id=paste0("met-d-", nmr_metabolites$GWAS_id)
nmr_metabolites_EG=nmr_metabolites[which(nmr_metabolites$Include=="yes"),]
nmr_metabolites_EG=nmr_metabolites_EG[-which(nmr_metabolites_EG$Group=="Glycolysis related metabolites"),]
mvmr_instruments=mv_extract_exposures(nmr_metabolites_EG$GWAS_id,pval_threshold=5e-8,clump_r2=0.01)
EG_snps=pull(mvmr_instruments, SNP)
expdat=extract_outcome_data(snps = EG_snps, outcome = nmr_metabolites_EG$GWAS_id)
names(expdat)=gsub("outcome", "exposure", names(expdat))
outcome_dat=read_outcome_data(snps=EG_snps, filename="Maternal_Effect_European_meta_NG2019.txt",snp_col="RSID", beta_col="beta", effect_allele_col="ea",
                               other_allele_col="nea", eaf_col="eaf", pval_col="p")
mvmr_harm=harmonise_data(expdat, outcome_dat, action=2) 
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
  effect_alleles[i,]=mvmr_harm[i,c("SNP", "effect_allele.exposure", "other_allele.exposure")]
}
merged=merge(effect_alleles, test, by="SNP")
reshaped_mvmr_data=merged
rs=reshaped_mvmr_data$SNP
bw_beta=outcome_dat$beta.outcome
bw_se=outcome_dat$se.outcome
find_beta_cols=grep("beta.exposure",names(reshaped_mvmr_data))
betaX_2=reshaped_mvmr_data[,find_beta_cols]
rf_2=colnames(betaX_2)
rs_2=reshaped_mvmr_data[, 1]
betaX_2.2=sapply(betaX_2, as.numeric )        
betaX_ivw_2=betaX_2.2/bw_se
bw_beta_ivw_2=bw_beta/bw_se
bw_nmr_input_2=new(Class = "mvMRInput", betaX = as.matrix(betaX_ivw_2), betaY=as.matrix(bw_beta_ivw_2), snps=rs_2, exposure=rf_2, outcome = "Birthweight")
ao=available_outcomes()
ao=ao[grep("met-d", ao$id),]
BMA_output=summarymvMR_SSS(bw_nmr_input_2, kmin=1, kmax=12, prior_prob=0.1, max_iter=100000)
mr.bw_BMA.out.step1=sss.report.mr.bma(BMA_output, top=10)
best.model.out.step1=sss.report.best.model(BMA_output, prior_sigma=0.5, top=10)

rename_traits=as.data.frame(best.model.out.step1)
rename_traits$`rf combination`=gsub("beta.exposure.","",rename_traits$`rf combination`)
rename_traits_mr=as.data.frame(mr.bw_BMA.out.step1)
rename_traits_mr$rf=gsub("beta.exposure.","",rename_traits_mr$rf)
ao=available_outcomes()
ao=ao[grep("met-d", ao$id),]
for (j in 1:nrow(ao))
{
  rename_traits$`rf combination`=gsub(ao$id[j], ao$trait[j], rename_traits$`rf combination`)
  rename_traits_mr$rf=gsub(ao$id[j], ao$trait[j], rename_traits_mr$rf)
}
best.model.out.updated=rename_traits
mr.bw_BMA.out.updated=rename_traits_mr
best.model.out.updated=apply(best.model.out.updated, 2, as.character)
mr.bw_BMA.out.updated=apply(mr.bw_BMA.out.updated, 2, as.character)

rm(list=ls())

## Kett
setwd()
nmr_metabolites=read_excel("nmr_metabolites.xlsx")
nmr_metabolites=nmr_metabolites[which(nmr_metabolites$Include=="yes"),]
colnames(nmr_metabolites)[1]="GWAS_id"
nmr_metabolites$GWAS_id=paste0("met-d-", nmr_metabolites$GWAS_id)
nmr_metabolites_EG=nmr_metabolites[which(nmr_metabolites$Include=="yes"),]
nmr_metabolites_EG=nmr_metabolites_EG[-which(nmr_metabolites_EG$Group=="Glycolysis related metabolites"),]
load("list_of_datasets")
EG_ids=list_of_datasets$id[-which(list_of_datasets$id%in%c("met-c-859", "met-c-920", "met-c-894", "met-c-849"))]
mvmr_instruments=mv_extract_exposures(EG_ids,pval_threshold=5e-8,clump_r2=0.01)
EG_snps=pull(mvmr_instruments, SNP)
expdat=extract_outcome_data(snps = EG_snps, outcome = nmr_metabolites_EG$GWAS_id)
names(expdat)=gsub("outcome", "exposure", names(expdat))
outcome_dat=read_outcome_data(snps=EG_snps, filename="Maternal_Effect_European_meta_NG2019.txt", snp_col="RSID", beta_col="beta", effect_allele_col="ea",
                               other_allele_col="nea", eaf_col="eaf", pval_col="p")
mvmr_harm=harmonise_data(expdat, outcome_dat, action=2) 

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
  effect_alleles[i,]=mvmr_harm[i,c("SNP", "effect_allele.exposure", "other_allele.exposure")]
}
merged=merge(effect_alleles, test, by="SNP")
outcome_dat=read_outcome_data(snps=merged$SNP,filename="Maternal_Effect_European_meta_NG2019.txt",snp_col="RSID", beta_col = "beta",effect_allele_col = "ea",
                              other_allele_col = "nea" , eaf_col = "eaf", pval_col = "p", phenotype_col = "Birthweight")
reshaped_mvmr_data=merged
rs=reshaped_mvmr_data$SNP
bw_beta=outcome_dat$beta.outcome
bw_se=outcome_dat$se.outcome
find_beta_cols=grep("beta.exposure",names(reshaped_mvmr_data))
betaX_2=reshaped_mvmr_data[,find_beta_cols]
rf_2=colnames(betaX_2)
rs_2=reshaped_mvmr_data[, 1]
betaX_2.2=sapply(betaX_2, as.numeric )        
betaX_ivw_2=betaX_2.2/bw_se
bw_beta_ivw_2=bw_beta/bw_se
bw_nmr_input_2=new(Class = "mvMRInput", betaX = as.matrix(betaX_ivw_2), betaY=as.matrix(bw_beta_ivw_2), snps=rs_2, exposure=rf_2, outcome = "Birthweight")
BMA_output=summarymvMR_SSS(bw_nmr_input_2, kmin=1, kmax=12, prior_prob=0.1, max_iter=100000)
mr.bw_BMA.out.step1=sss.report.mr.bma(BMA_output, top=10)
best.model.out.step1=sss.report.best.model(BMA_output, prior_sigma=0.5, top=10)

rename_traits=as.data.frame(best.model.out.step1)
rename_traits$`rf combination`=gsub("beta.exposure.","",rename_traits$`rf combination`)
rename_traits_mr=as.data.frame(mr.bw_BMA.out.step1)
rename_traits_mr$rf=gsub("beta.exposure.","",rename_traits_mr$rf)
ao=available_outcomes()
ao=ao[grep("met-d", ao$id),]
for (j in 1:nrow(ao))
{
  rename_traits$`rf combination`=gsub(ao$id[j], ao$trait[j], rename_traits$`rf combination`)
  rename_traits_mr$rf=gsub(ao$id[j], ao$trait[j], rename_traits_mr$rf)
}
best.model.out.updated=rename_traits
mr.bw_BMA.out.updated=rename_traits_mr
best.model.out.updated=apply(best.model.out.updated, 2, as.character)
mr.bw_BMA.out.updated=apply(mr.bw_BMA.out.updated, 2, as.character)

rm(list=ls())

## Multi
ket_gwas=read_csv("ket_gwas")
ket_gwas=ket_gwas[-which(ket_gwas$class%in%c("Glycolysis related metabolites")),]
ao=available_outcomes()
ao_1=ao[which(ao$year==2016&ao$author=="Kettunen"&ao$category=="Metabolites"&ao$population=="European"),]
ket_gwas=ao_1[which(ao_1$trait%in%ket_gwas$id.exposure),]
# snps from other GWAS
other=c("ieu-a-300","ieu-a-299", "ieu-a-302", "ebi-a-GCST005186", "ebi-a-GCST004939", "ieu-a-1105")
ldl_c=ao[grep("ieu-a-300", ao$id),]
hdl_c=ao[grep("ieu-a-299", ao$id),]
TG=ao[grep("ieu-a-302", ao$id),]
Serum_creatinine=ao[grep("ieu-a-1105", ao$id),]
find_snps=rbind(ket_gwas, ldl_c, hdl_c, TG, Serum_creatinine)
multi_snps=mv_extract_exposures(find_snps$id, pval_threshold=5e-8,clump_r2=0.01)
EG_snps=pull(multi_snps, SNP)
expdat=extract_outcome_data(snps = EG_snps, outcome = nmr_metabolites_EG$GWAS_id)
names(expdat)=gsub("outcome", "exposure", names(expdat))
outcome_dat=read_outcome_data(snps=EG_snps, filename="UKBB_birthweight",snp_col="RSID", beta_col="beta", effect_allele_col="ea",
                               other_allele_col="nea", eaf_col="eaf", pval_col="p")
mvmr_harm=harmonise_data(expdat, outcome_dat, action=2) 
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
  effect_alleles[i,]=mvmr_harm[i,c("SNP", "effect_allele.exposure", "other_allele.exposure")]
}
merged=merge(effect_alleles, test, by="SNP")
outcome_dat=read_outcome_data(snps=merged$SNP,filename="UKBB_birthweight",snp_col="RSID", beta_col = "beta",effect_allele_col = "ea",
                              other_allele_col = "nea" , eaf_col = "eaf", pval_col = "p", phenotype_col = "Birthweight")
reshaped_mvmr_data=merged
rs=reshaped_mvmr_data$SNP
bw_beta=outcome_dat$beta.outcome
bw_se=outcome_dat$se.outcome
find_beta_cols=grep("beta.exposure",names(reshaped_mvmr_data))
betaX_2=reshaped_mvmr_data[,find_beta_cols]
rf_2=colnames(betaX_2)
rs_2=reshaped_mvmr_data[, 1]
betaX_2.2=sapply(betaX_2, as.numeric )        
betaX_ivw_2=betaX_2.2/bw_se
bw_beta_ivw_2=bw_beta/bw_se
bw_nmr_input_2=new(Class = "mvMRInput", betaX = as.matrix(betaX_ivw_2), betaY=as.matrix(bw_beta_ivw_2), snps=rs_2, exposure=rf_2, outcome = "Birthweight")
BMA_output=summarymvMR_SSS(bw_nmr_input_2, kmin=1, kmax=12, prior_prob=0.1, max_iter=100000)
mr.bw_BMA.out.step1=sss.report.mr.bma(BMA_output, top=10)
best.model.out.step1=sss.report.best.model(BMA_output, prior_sigma=0.5, top=10)

rename_traits=as.data.frame(best.model.out.step1)
rename_traits$`rf combination`=gsub("beta.exposure.","",rename_traits$`rf combination`)
rename_traits_mr=as.data.frame(mr.bw_BMA.out.step1)
rename_traits_mr$rf=gsub("beta.exposure.","",rename_traits_mr$rf)
ao=available_outcomes()
ao=ao[grep("met-d", ao$id),]
for (j in 1:nrow(ao))
{
  rename_traits$`rf combination`=gsub(ao$id[j], ao$trait[j], rename_traits$`rf combination`)
  rename_traits_mr$rf=gsub(ao$id[j], ao$trait[j], rename_traits_mr$rf)
}
best.model.out.updated=rename_traits
mr.bw_BMA.out.updated=rename_traits_mr
best.model.out.updated=apply(best.model.out.updated, 2, as.character)
mr.bw_BMA.out.updated=apply(mr.bw_BMA.out.updated, 2, as.character)

rm(list=ls())
