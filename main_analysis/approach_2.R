library(readr); library(ggplot2); library(data.table);library(TwoSampleMR);library(stringr);library(readxl)
source("VZ_summary_mvMR_SSS_function.R")
source("VZ_summary_mvMR_BF_function.R")

# II. Kettunen SNPs only

# 1. Select genetic instruments, p-value<5x10^-8, Rsq<0.01
nmr_metabolites=read_excel("nmr_metabolites.xlsx")
nmr_metabolites=nmr_metabolites[which(nmr_metabolites$Include=="yes"),]
colnames(nmr_metabolites)[1]="GWAS_id"
nmr_metabolites$GWAS_id=paste0("met-d-", nmr_metabolites$GWAS_id)
ao=available_outcomes()
ao=ao[which(ao$year==2016&ao$author=="Kettunen"&ao$category=="Metabolites"&ao$population=="European"),]
ao_1=ao[which(ao$trait%in%nmr_metabolites$`Biomarker name`),]
name_mismatch=c("18:2, linoleic acid (LA)", "Free cholesterol", "3-hydroxybutyrate", "Apolipoprotein A-I", "22:6, docosahexaenoic acid", "Mono-unsaturated fatty acids",
                "Omega-7, omega-9 and saturated fatty acids", "Phosphatidylcholine and other cholines", "Total lipids in chylomicrons and largest VLDL particles", "Serum total triglycerides")
ao_1=rbind(ao_1, ao[which(ao$trait%in%name_mismatch),] )
list_of_datasets=as.data.frame(ao_1)

# Kettunen dataset metabolites missing: Acetone, Glycine, Total esterified cholesterol, Total phospholipids in lipoprotein particles

ket_ids=as.character(list_of_datasets$id)
uvmr_instruments=extract_instruments(ket_ids,p1=5e-8,clump=T,r2=0.01)


# 2. Extract SNP-exposure data for select SNPs
exposure_data=uvmr_instruments

# 3. Extract snp-outcome data for select SNPs, univariate instruments

# load outcome data - Warrington GWAS of maternal genetic effects on birthweight, adjusted for fetal genotype
outcome_dat=read_outcome_data(snps=exposure_data$SNP, filename="Maternal_Effect_European_meta_NG2019.txt",snp_col="RSID", beta_col="beta",
                              effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p")

 # 4. Harmonise snp-exposure and snp-outcome data
uvmr_harm=harmonise_data(exposure_data, outcome_dat, action=2)

# 5. Run Wald ratio or IVW
estimates=mr(uvmr_harm, method_list=c("mr_wald_ratio", "mr_ivw"))
colnames(estimates)[1]="id"
ao=available_outcomes()
ao_1=ao[which(ao$year==2016&ao$author=="Kettunen"&ao$category=="Metabolites"&ao$population=="European"),]
ao_1=ao_1[which(ao_1$id%in%estimates$id),c("id", "trait")]
estimates=merge(estimates, ao_1, by="id")

# 6. Run BMA-MVMR on Bluepebble
mvmr_instruments=mv_extract_exposures(unlist(list_of_datasets$id),pval_threshold=5e-8,clump_r2=0.01)
ket_snps=pull(mvmr_instruments, SNP)
expdat=extract_outcome_data(snps = ket_snps, outcome = nmr_metabolites$GWAS_id)
names(expdat)=gsub("outcome", "exposure", names(expdat))
outcome_dat=read_outcome_data(snps=ket_snps, filename="Maternal_Effect_European_meta_NG2019.txt",
                              snp_col="RSID", beta_col="beta", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p")
mvmr_harm=harmonise_data(expdat, outcome_dat, action=2) 

setDT(mvmr_harm)
test=reshape(mvmr_harm, timevar="id.exposure", idvar=c("SNP"), direction="wide")
test=as.matrix(test)
col=as.numeric(grep("exposure|SNP", colnames(test)))
test=test[,col]
length(unique(test[,1]))
snps=unique(mvmr_harm$SNP)
effect_alleles=as.data.frame(matrix(1:length(snps),nrow=length(snps),ncol=3))
colnames(effect_alleles)=c("SNP", "eaf.outcome", "remove")
for (i in 1:length(snps))
{
  effect_alleles[i,]=mvmr_harm[i,c("SNP", "eaf.outcome", "remove")]
}
merged=merge(effect_alleles, test, by="SNP")
outcome_dat=read_outcome_data(snps=merged$SNP,filename="Maternal_Effect_European_meta_NG2019.txt", snp_col="RSID", beta_col = "beta",
                              effect_allele_col = "ea", other_allele_col = "nea" , eaf_col = "eaf", pval_col = "p")
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
qchisq(0.05, 111, lower.tail=F)

# remove: rs34346326, rs41272659, rs12325419, rs118182637, rs2079742, rs4410345, rs4149307
# rs10211524, rs2657879, rs6729711, rs145717049

# Q-statistics, df~= number of SNPs
model_index=best.model.out.step1[1,1]
title=ao$trait[which(ao$id==gsub("beta.exposure.","",model_index))]
Q = matrix(ncol=1, nrow=length(bw_beta_ivw_2))
betaX_model=as.matrix(betaX_ivw_2[,unlist(model_index, ",")])
sigma_vec=rep(0.5, ncol(betaX_model))
H_fm=betaX_model%*% solve(t(betaX_model) %*% betaX_model + sigma_vec^{-2} ) %*% t(betaX_model)
predicted_bw=H_fm %*% bw_beta_ivw_2
Q=(bw_beta_ivw_2-predicted_bw)^2
maxQ=apply(Q, MARGIN=1, FUN=max)
sort.ix = sort.int(maxQ, decreasing=TRUE, index.return=TRUE)
Q_tab=cbind(rs_2,Q) 
qchisq(0.05, 111, lower.tail=F)

write.csv(Q_tab, "kett_all_qstats.csv")
# remove: rs34346326, rs41272659, rs12325419, rs118182637, rs2079742, rs4410345, rs4149307
# rs10211524, rs2657879, rs6729711, rs145717049

png("q_plots_kett.png", height = 1000, width=1400)
df = data.frame(x=predicted_bw, y =bw_beta_ivw_2, Q = Q, rs=rs_2[which(rs_2%in%outcome_dat$SNP)])
p = ggplot(df, aes(x , y)) +  geom_point(aes(colour = Q), size =4) + scale_colour_gradientn(colours =c("#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"))  +
  labs(x = "Predicted beta bwt", y="Observed beta bwt", colour="Q") +
  geom_hline(yintercept = 0, linetype="dotted") + geom_vline(xintercept = 0, linetype="dotted") +
  geom_text(aes(label=ifelse(Q>6,as.character(rs[which(rs%in%outcome_dat$SNP)]),'')),position=position_jitter(width=.35,height=0.4),hjust=0.5, vjust=-1, size=5) +
  theme(axis.text.x = element_text(size = 4), axis.text.y = element_text(size = 4), axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 16), legend.text=element_text(size=16),legend.title=element_text(size=18)) +
  ggtitle(title)
p
dev.off()

# Cooks D
sigma_vec=rep(0.5, ncol(betaX_model))
cD=as.numeric(cooksD(bw_beta_ivw_2, betaX_model,sigma_vec)$cooksD)
cD_thresh=cooksD(bw_beta_ivw_2, betaX_model, sigma_vec)$cooksD_thresh
cooksD_tab=data.frame(rs_2,cD)

png("CD_plots_kett.png", height = 1000, width=1400)
df=data.frame(x=predicted_bw, y =bw_beta_ivw_2, cD = cD, rs=rs_2[which(rs_2%in%outcome_dat$SNP)])
p=ggplot(df, aes(x, y)) + geom_point(aes(colour = cD), size =4) +
  scale_colour_gradientn(colours = c("white", "orange", "red", "darkred"))  +
  labs(x = "Predicted beta bwt", y="Observed beta bwt", colour="Cooks D") + geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted")  +
  geom_text(aes(label=ifelse(cD>cD_thresh,as.character(rs[which(rs%in%outcome_dat$SNP)]),'')),hjust=0.5, vjust=-1, size=5)+
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13), axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18), legend.text=element_text(size=16),legend.title=element_text(size=18))+
  ggtitle(title)
p
dev.off()
