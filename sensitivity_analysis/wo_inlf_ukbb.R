# Excluding influential variants, using q-statistic
setwd()
library(readr); library(ggplot2); library(data.table);library(TwoSampleMR);library(stringr);library(readxl)
library(tidyverse); library(ggforestplot); library(mr.raps)
source("VZ_summary_mvMR_SSS_function.R")
source("VZ_summary_mvMR_BF_function.R")
`%!in%`=Negate(`%in%`)

mvmr_harm_wo_infl=read_csv("mvmr_harm.csv")
infl=c("rs11601507", "rs35004366", "rs560887", "rs4844599", "rs12369443", "rs10008637", "rs2467853", "rs1958029", "rs61705884",
       "rs1498694", "rs34362588", "rs534417", "rs1274961", "rs3768321", "rs35823804", "rs40270", "rs4656292", "rs17369578", "rs12524288",
       "rs143699220", "rs7150045", "rs11976955", "rs4813543", "rs9389268", "rs218674", "rs77960347", "rs10062079", "rs34707604", "rs56139160",
       "rs145679432", "rs7078003", "rs1965869", "rs11564722", "rs12206654", "rs2631367", "rs113312468")
mvmr_harm_wo_infl=mvmr_harm_wo_infl[-which(mvmr_harm_wo_infl$SNP%in%infl),]
setDT(mvmr_harm_wo_infl)
test=reshape(mvmr_harm_wo_infl, timevar="id.exposure", idvar=c("SNP"), direction="wide")
test=as.matrix(test)
col=as.numeric(grep("exposure|SNP", colnames(test)))
test=test[,col]
length(unique(test[,1]))
snps=unique(mvmr_harm_wo_infl$SNP)
effect_alleles=as.data.frame(matrix(1:length(snps),nrow=length(snps),ncol=3))
colnames(effect_alleles)=names(mvmr_harm_wo_infl)[c(1,9,10)]
for (i in 1:length(snps))
{
  effect_alleles[i,1:3]=mvmr_harm_wo_infl[i,c(1,9,10)]
}
merged=merge(effect_alleles, test, by="SNP")
outcome_dat=read_outcome_data(snps=merged$SNP,filename="UKBB_birthweight",snp_col="RSID", beta_col = "beta",effect_allele_col = "ea",
                              other_allele_col = "nea" , eaf_col = "eaf", pval_col = "p",samplesize_col ="n_ownBW",phenotype_col = "Birthweight")

test=reshape(mvmr_harm_wo_infl, timevar="id.exposure", idvar=c("SNP"), direction="wide")
test=as.matrix(test)
col=as.numeric(grep("exposure|SNP", colnames(test)))
test=test[,col]
length(unique(test[,1]))
snps=unique(mvmr_harm_wo_infl$SNP)
effect_alleles=as.data.frame(matrix(1:length(snps),nrow=length(snps),ncol=3))
colnames(effect_alleles)=names(mvmr_harm_wo_infl)[c(1,9,10)]
for (i in 1:length(snps))
{
  effect_alleles[i,1:3]=mvmr_harm_wo_infl[i,c(1,9,10)]
}
merged=merge(effect_alleles, test, by="SNP")

bw_beta=outcome_dat$beta.outcome
bw_se=outcome_dat$se.outcome
find_beta_cols=grep("beta.exposure",names(merged))
betaX_2=merged[,find_beta_cols]
rf_2=colnames(betaX_2)
rs_2=merged[, 1]
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
