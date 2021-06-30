# Excluding influential variants, using q-statistic

setwd()
devtools::install_github("NightingaleHealth/ggforestplot")
library(readr); library(ggplot2); library(data.table);library(TwoSampleMR);library(stringr);library(readxl)
library(tidyverse); library(ggforestplot); library(mr.raps)
source("VZ_summary_mvMR_SSS_function.R")
source("VZ_summary_mvMR_BF_function.R")
`%!in%`=Negate(`%in%`)

mvmr_harm_wo_infl=read_csv("kett_mvmr_harm.csv")
influ=c("rs34346326", "rs41272659", "rs12325419", "rs118182637", "rs2079742", "rs4410345", "rs4149307",
        "rs10211524", "rs2657879", "rs6729711", "rs145717049")
mvmr_harm_wo_infl=mvmr_harm_wo_infl[-which(mvmr_harm_wo_infl$SNP%in%influ),]
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
