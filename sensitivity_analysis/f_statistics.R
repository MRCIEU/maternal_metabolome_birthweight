## F statistics
# Univariate analysis -- mean F-statistic
# Multivariate analysis -- conditional F-statistic

library(readr);library(TwoSampleMR);library(data.table);library(readxl);library(remotes);library(stringr)
ao=available_outcomes()
`%!in%`=Negate(`%in%`)
install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE, force=T)

pheno_cov=read.delim("pheno_cov")
nmr_metabolites_UKBB=read_xlsx("nmr_metabolites")
nmr_metabolites_UKBB=nmr_metabolites_UKBB[which(nmr_metabolites_UKBB$Include=="yes"),]
colnames(ao)[2]="Biomarker name"
ao1=ao[which(ao$author=="Borges CM"&ao$`Biomarker name`%in%nmr_metabolites_UKBB$`Biomarker name`),]
nmr_metabolites_UKBB=merge(nmr_metabolites_UKBB, ao1, by="Biomarker name")
nmr_metabolites_UKBB$id=sub("met-d-", "", nmr_metabolites_UKBB$id)
pheno_cov=pheno_cov[which(rownames(pheno_cov)%in%nmr_metabolites_UKBB$id), which(colnames(pheno_cov)%in%nmr_metabolites_UKBB$id)]

## UKBB
# mean f-statistic
models=read_csv("models")
uvmr_harm=read_csv("uvmr_harm")
uvmr_harm$pre_fstat=(uvmr_harm$beta.exposure^2)/(uvmr_harm$se.exposure^2)
mean_fstat_ukbb=data.frame(NA,NA)
colnames(mean_fstat_ukbb)=c("model", "f-stat")
for (i in models$gwas_id)
{
  j=which(models$gwas_id==i)
  mean_fstat_ukbb[j,]=c(paste(i), mean(uvmr_harm$pre_fstat[which(uvmr_harm$id.exposure==i)]))
}

rm(list=ls())

## Kettunen
# mean f-statistic
models=read_csv("models")
uvmr_harm=read_csv("uvmr_harm")
uvmr_harm$pre_fstat=(uvmr_harm$beta.exposure^2)/(uvmr_harm$se.exposure^2)
mean_fstat_kett=data.frame(NA,NA)
colnames(mean_fstat_kett)=c("model", "f-stat")
for (i in na.omit(models)$gwas_id)
{
  j=which(na.omit(models)$gwas_id==i)
  mean_fstat_kett[j,]=c(paste(i), mean(uvmr_harm$pre_fstat[which(uvmr_harm$id.exposure==i)]))
}

# conditional f-statistic for multiple
combinations=str_split(models$`rf combination`[grep(",", models$`rf combination`)], ",")
ao1=ao[which(ao$author=="Borges CM"&ao$trait%in%unlist(str_split(models$`rf combination`[grep(",", models$`rf combination`)], ","))),]

list_c=list()
for (i in 1:length(combinations))
{
  a=unlist(combinations[i])
  b=ao1$id[which(ao1$trait%in%a)]
  c=gsub("met-d-","",b)
  mvmr_instr=mv_extract_exposures(b, clump_r2=0.01)
  outcome_dat=read_outcome_data(snps=mvmr_instr$SNP, filename="UKBB_birthweight",snp_col="RSID", beta_col="beta", effect_allele_col="ea",
                                 other_allele_col="nea", eaf_col="eaf", pval_col="p", samplesize_col="n_ownBW")
  mvmr_harm=mv_harmonise_data(mvmr_instr, outcome_dat, harmonise_strictness=2)
  pheno_cov_tmp=as.matrix(pheno_cov[which(rownames(pheno_cov)%in%c), which(colnames(pheno_cov)%in%c)])
  nmr_tmp=nmr_metabolites_UKBB[which(nmr_metabolites_UKBB$`Biomarker name`%in%combinations[[i]]),]
  
  BGX=mvmr_harm$exposure_beta
  seBGX=mvmr_harm$exposure_se

  BGY=mvmr_harm$outcome_beta
  seBGY=mvmr_harm$outcome_se

  correlation=MVMR::phenocov_mvmr(pheno_cov_tmp, seBGX)
  formatted=MVMR::format_mvmr(BXGs=BGX, seBXGs=seBGX, seBYG=seBGY, BYG=BGY, RSID=rownames(mvmr_harm$exposure_beta))
  list_c[[i]]=MVMR::strength_mvmr(r_input=formatted, correlation)
}

combi=unlist(combinations, recursive = F)
combi=data.frame(combi[seq(1,13, by=2)], combi[seq(2,14, by=2)], (unlist(list_c)[seq(1, 13, by=2)]), (unlist(list_c)[seq(2, 14, by=2)]))
colnames(combi)=c("exposure_1", "exposure_2", "F_stat_1", "F_stat_2")

rm(list=ls())

## Multi
# mean f-statistic
models=read_csv("models")

# note - some metabolites suggested only available in UKBB
exposure_data=extract_instruments(models$gwas_id, p1=5e-8, clump=T, r2=0.01)
outcome_dat=read_outcome_data(snps=exposure_data$SNP, filename="UKBB_birthweight",snp_col="RSID", beta_col="beta", effect_allele_col="ea",
                              other_allele_col="nea", eaf_col="eaf", pval_col="p", samplesize_col="n_ownBW")
uvmr_harm=harmonise_data(exposure_data, outcome_dat, action=2)
uvmr_harm$pre_fstat=(uvmr_harm$beta.exposure^2)/(uvmr_harm$se.exposure^2)
mean_fstat_multi=data.frame(NA,NA)
colnames(mean_fstat_multi)=c("model", "f-stat")
for (i in models$gwas_id)
{
  j=which(models$gwas_id==i)
  mean_fstat_multi[j,]=c(paste(i), mean(uvmr_harm$pre_fstat[which(uvmr_harm$id.exposure==i)]))
}

