library(dplyr);library(readr);library(TwoSampleMR);library(stringr)
`%!in%`=Negate(`%in%`)

## UKBB
# UVMR
setwd()
exposure_data=read_csv("uvmr_instruments_ukbb.csv")[,-1]
outcome_dat=read_outcome_data(snps=uvmr_instruments$SNP, filename="Maternal_Effect_European_meta_NG2019.txt",
                              snp_col="RSID", beta_col="beta", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p")
uvmr_harm=harmonise_data(exposure_data, outcome_dat, action=2)
estimates=mr(uvmr_harm, method_list=c("mr_wald_ratio", "mr_ivw"))
ao=available_outcomes()
ao_1=ao[which(ao$id%in%estimates$id.exposure),]
for (j in 1:nrow(ao_1))
{
  estimates$id.exposure=gsub(ao_1$id[j], ao_1$trait[j], estimates$id.exposure)
}
models=as.data.frame(read_csv("ukbb_best_model_out_default.csv")[,-1])
estimates=estimates[order(estimates$id.exposure),]
models=models[order(models$`rf combination`),]
models[,4:5]=estimates[which(estimates$id.exposure%in%models$`rf combination`),c("b","se")]
#models[,6]=do.call(paste, c("(",round(models[,4]-1.96*models[,5],3),",",round(models[,4]+1.96*models[,5],3),")", sep = ""))
models[,6]=round(models[,4]-1.96*models[,5],3)
models[,7]=round(models[,4]+1.96*models[,5],3)
models$b=round(models$b, 3)
models$se=round(models$se, 3)


## Kettunen
# UVMR
setwd("")
estimates=read.csv("uvmr_estimates.csv")
models=read_csv("kett_best_model_out_default_2.csv")[,-1]
uvmr_subset=estimates[which(estimates$trait%in%models$`rf combination`[which(1:10%!in%grep(",", models$`rf combination`))]),c("trait","b","se")]
uvmr_subset[,4]=round(uvmr_subset[,2]-1.96*uvmr_subset[,3],3)
uvmr_subset[,5]=round(uvmr_subset[,2]+1.96*uvmr_subset[,3], 3)
uvmr_subset[,6]=paste0("(",uvmr_subset[,4],", ", uvmr_subset[,5],")")

# MVMR
ao1=ao[which(ao$author=="Borges CM"&ao$trait%in%unlist(str_split(models$`rf combination`[grep(",", models$`rf combination`)], ","))),]
combinations=str_split(models$`rf combination`[grep(",", models$`rf combination`)], ",")
list_c=list()
for (i in 1:length(combinations))
{
  a=unlist(combinations[i])
  b=ao1$id[which(ao1$trait%in%a)]
  mvmr_instr=mv_extract_exposures(b, clump_r2=0.01)
  snps=pull(mvmr_instr, SNP)
  outcome_dat=read_outcome_data(snps=snps, filename="Maternal_Effect_European_meta_NG2019.txt",
                                snp_col="RSID", beta_col="beta", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p")
  mvmr_harm=mv_harmonise_data(mvmr_instr, outcome_dat, harmonise_strictness=2)
  list_c[i]=mv_basic(mvmr_harm, pval_threshold = 5e-08)
}

# create subsets of combinations, extract instruments, extract outcome, harmonise and calculate estimate
for (i in 1:length(list_c))
{
  print(list_c[[i]]$id.exposure)
  print(round(list_c[[i]]$b,3))
  print(paste0("(",round(list_c[[i]]$b-1.96*list_c[[i]]$se,3),", ",round(list_c[[i]]$b+1.96*list_c[[i]]$se,3),")"))
  print(round(list_c[[i]]$se,3))
}

