devtools::install_github("NightingaleHealth/ggforestplot")
library(readr); library(ggplot2); library(data.table);library(TwoSampleMR);library(stringr);library(readxl)
library(tidyverse); library(ggforestplot)
source("VZ_summary_mvMR_SSS_function.R")
source("VZ_summary_mvMR_BF_function.R")
`%!in%`=Negate(`%in%`)

# I. UKBB SNPs only

# 1. Select genetic instruments, p-value<5x10^-8, Rsq<0.01
nmr_metabolites_UKBB=read_excel("nmr_metabolites")
ao=available_outcomes()
trait_linkage=ao[grep("met-d", ao$id),]
nmr_metabolites_UKBB=nmr_metabolites_UKBB[which(nmr_metabolites_UKBB$Include=="yes"),]
nmr_metabolites_UKBB$Name=paste0("met-d-", nmr_metabolites_UKBB$Name)
names(nmr_metabolites_UKBB)[which(names(nmr_metabolites_UKBB)=="Name")]="GWAS_id"

# 2. Extract SNP-exposure data for select SNPs, univariate instruments
exposure_data=extract_instruments(nmr_metabolites_UKBB$GWAS_id, p1=5e-8, clump=T, r2=0.01)

# 3. Extract snp-outcome data for select SNPs, univariate instruments
# load outcome data - Warrington GWAS of maternal genetic effects on birthweight, adjusted for fetal genotype
outcome_dat=read_outcome_data(snps=exposure_data$SNP, filename="UKBB_birthweight",
                              snp_col="RSID", beta_col="beta", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p", samplesize_col="n_ownBW")

# 4. Harmonise snp-exposure and snp-outcome data
uvmr_harm=harmonise_data(exposure_data, outcome_dat, action=2)

# 5. Run Wald ratio or IVW
estimates=mr(uvmr_harm, method_list=c("mr_wald_ratio", "mr_ivw"))
ao_1=ao[which(ao$id%in%estimates$id.exposure),]
for (j in 1:nrow(ao_1))
{
  estimates$id.exposure=gsub(ao_1$id[j], ao_1$trait[j], estimates$id.exposure)
}
ao_1=ao[which(ao$year=="2020"&ao$author=="Borges CM"&ao$population=="European"),]
colnames(estimates)[1]="trait"
estimates=merge(estimates, ao_1, by="trait")


# 6. Run BMA-MVMR
# multivariate instruments
mvmr_instruments_ukbb=mv_extract_exposures(nmr_metabolites_UKBB$GWAS_id, clump_r2=0.01)

nmr_metabolites_UKBB=read_excel("nmr_metabolites")
setDT(mvmr_instruments_ukbb)
data_wide=dcast(mvmr_instruments_ukbb,SNP~id.exposure,value.var=list(names(mvmr_instruments_ukbb)[6:9]))
head(mvmr_instruments_ukbb)

# make a matrix of snp, effect allele and other allele to merge
snps=unique(mvmr_instruments_ukbb$SNP)
effect_alleles=as.data.frame(matrix(1:length(snps),nrow=length(snps),ncol=3))
colnames(effect_alleles)=names(mvmr_instruments_ukbb)[c(1,4,5)]
for (i in 1:length(snps))
{
  effect_alleles[i,1:3]=mvmr_instruments_ukbb[i,c(1,4,5)]
}
merged=merge(effect_alleles, data_wide, by="SNP")
head(merged)
reshaped_mvmr_data=merged
outcome_dat=read_outcome_data(snps=mvmr_instruments_ukbb$SNP,filename="UKBB_birthweight",snp_col="RSID", beta_col = "beta",
                              effect_allele_col = "ea",other_allele_col = "nea" , eaf_col = "eaf", pval_col = "p",samplesize_col ="n_ownBW",phenotype_col = "Birthweight")
mvmr_harm=harmonise_data(mvmr_instruments_ukbb, outcome_dat, action=2) 

setDT(mvmr_harm)
test=reshape(mvmr_harm, timevar="id.exposure", idvar=c("SNP"), direction="wide")
test=as.matrix(test)
col=as.numeric(grep("exposure|SNP", colnames(test)))
test=test[,col]
length(unique(test[,1]))
snps=unique(mvmr_harm$SNP)
effect_alleles=as.data.frame(matrix(1:length(snps),nrow=length(snps),ncol=3))
colnames(effect_alleles)=names(mvmr_harm)[c(1,9,10)]
for (i in 1:length(snps))
{
  effect_alleles[i,1:3]=mvmr_harm[i,c(1,9,10)]
}
merged=merge(effect_alleles, test, by="SNP")
outcome_dat=read_outcome_data(snps=merged$SNP,filename="UKBB_birthweight",snp_col="RSID", beta_col = "beta",
                              effect_allele_col = "ea",other_allele_col = "nea" , eaf_col = "eaf", pval_col = "p",samplesize_col ="n_ownBW",phenotype_col = "Birthweight")
# to add in trait names for biological interpretation
nmr_metabolites_UKBB=nmr_metabolites_UKBB[which(nmr_metabolites_UKBB$Include=="yes"),]
nmr_metabolites_UKBB$Name=paste0("met-d-", nmr_metabolites_UKBB$Name)
names(nmr_metabolites_UKBB)[which(names(nmr_metabolites_UKBB)=="Biomarker name")]="GWAS_id"
ao=available_outcomes()
ao=ao[grep("met-d", ao$id),]
ao=ao[which(ao$id%in%nmr_metabolites_UKBB$Name),]

# create bma-mvmr input class
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

# default model prior prob = 0.1, prior sigma = 0.5
BMA_output=summarymvMR_SSS(bw_nmr_input_2, kmin=1, kmax=12, prior_prob=0.1, max_iter=100000)
mr.bw_BMA.out.step1=sss.report.mr.bma(BMA_output, top=10)
best.model.out.step1=sss.report.best.model(BMA_output, prior_sigma=0.5, top=10)

rename_traits=as.data.frame(best.model.out.step1)
rename_traits$`rf combination`=gsub("beta.exposure_","",rename_traits$`rf combination`)
rename_traits_mr=as.data.frame(mr.bw_BMA.out.step1)
rename_traits_mr$rf=gsub("beta.exposure_","",rename_traits$rf)
for (j in 1:nrow(ao))
{
  rename_traits$`rf combination`=gsub(ao$id[j], ao$trait[j], rename_traits$`rf combination`)
  rename_traits_mr$rf=gsub(ao$id[j], ao$trait[j], rename_traits_mr$rf)
}
best.model.out.updated=rename_traits
mr.bw_BMA.out.updated=rename_traits_mr
best.model.out.updated=apply(best.model.out.updated, 2, as.character)
mr.bw_BMA.out.updated=apply(mr.bw_BMA.out.updated, 2, as.character)

# Q-statistics, df~= number of SNPs
model_index=best.model.out.step1[1,1]
Q = matrix(ncol=1, nrow=length(bw_beta_ivw_2))
betaX_model=as.matrix(betaX_ivw_2[,unlist(model_index, ",")])
title=model_index
sigma_vec=rep(0.5, ncol(betaX_model))
H_fm=betaX_model%*% solve(t(betaX_model) %*% betaX_model + sigma_vec^{-2} ) %*% t(betaX_model)
predicted_bw=H_fm %*% bw_beta_ivw_2
Q=(bw_beta_ivw_2-predicted_bw)^2
maxQ=apply(Q, MARGIN=1, FUN=max)
sort.ix = sort.int(maxQ, decreasing=TRUE, index.return=TRUE)
Q_tab=cbind(rs_2,Q) 
qchisq(0.05, 610, lower.tail=F)

# remove: rs11601507, rs35004366, rs560887, rs4844599, rs12369443, rs10008637, rs2467853, rs1958029, rs61705884
# rs1498694, rs34362588, rs534417, rs1274961, rs3768321, rs35823804, rs40270, rs4656292, rs17369578, rs12524288
# rs143699220, rs7150045, rs11976955, rs4813543, rs9389268, rs218674, rs77960347, rs10062079, rs34707604, rs56139160
# rs145679432, rs7078003, rs1965869, rs11564722, rs12206654, rs2631367, rs113312468

