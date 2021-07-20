## additional plots
library(readr);library(data.table);library(TwoSampleMR)

## UKBB data
ukbb_mvmr_harm=read_csv("ukbb_mvmr_harm.csv")
setDT(ukbb_mvmr_harm)
test=reshape(ukbb_mvmr_harm, timevar="id.exposure", idvar=c("SNP"), direction="wide")
test=as.matrix(test)
col=as.numeric(grep("exposure|SNP", colnames(test)))
test=test[,col]
length(unique(test[,1]))
snps=unique(ukbb_mvmr_harm$SNP)
effect_alleles=as.data.frame(matrix(1:length(snps),nrow=length(snps),ncol=3))
colnames(effect_alleles)=c("SNP", "eaf.outcome", "remove")
for (i in 1:length(snps))
{
        effect_alleles[i,]=ukbb_mvmr_harm[i,c("SNP", "eaf.outcome", "remove")]
}
merged=merge(effect_alleles, test, by="SNP")
outcome_dat=read_outcome_data(snps=merged$SNP,filename="UKBB_birthweight", snp_col="RSID", beta_col = "beta",
                              effect_allele_col = "ea", other_allele_col = "nea" , eaf_col = "eaf", pval_col = "p", phenotype_col = "Birthweight")
beta_x_leu=as.numeric(merged$`beta.exposure.met-d-Leu`)
se_x_leu=as.numeric(merged$`se.exposure.met-d-Leu`)
beta_y_leu=as.numeric(outcome_dat$beta.outcome)
se_y_leu=as.numeric(outcome_dat$se.outcome)

beta_x_glu=as.numeric(merged$`beta.exposure.met-d-Glucose`)
se_x_glu=as.numeric(merged$`se.exposure.met-d-Glucose`)
beta_y_glu=as.numeric(outcome_dat$beta.outcome)
se_y_glu=as.numeric(outcome_dat$se.outcome)

png("snp-glucose+leucine", width = 800, height = 800)
plot(beta_x_leu, beta_y_leu, pch=16, ylab="SNP-BWT assoc.", xlab = "SNP-metabolite assoc.", cex.lab=1.2, cex.axis=1.5,
     cex=1.5, xlim=c(-0.15, 0.15), ylim=c(-0.1, 0.1))
arrows(beta_x_leu, beta_y_leu-(se_y_leu*1.96), beta_x_leu,
       beta_y_leu+(se_y_leu*1.96), length=0, angle=90, code=3, lty=2)
arrows(beta_x_leu-(se_x_leu*1.96), beta_y_leu, beta_x_leu+
               (se_x_leu*1.96), beta_y_leu, length=0, angle=90, code=3, lty=2)
points(beta_x_glu, beta_y_glu, col="red",pch=2)
arrows(beta_x_glu, beta_y_glu-(se_y_glu*1.96), beta_x_glu,
       beta_y_glu+(se_y_glu*1.96), length=0, angle=90, code=3, lty=2,col="red")
arrows(beta_x_glu-(se_x_glu*1.96), beta_y_glu, beta_x_glu+
               (se_x_glu*1.96), beta_y_glu, length=0, angle=90, code=3, lty=2,col="red")
legend("bottomleft", 
       legend = c("Leucine", "Glucose"), 
       col = c("black","red"), 
       pch = c(16,2), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))
dev.off()

## Kettunen data
kett_mvmr_harm=read_csv("kett_mvmr_harm")
setDT(kett_mvmr_harm)
test=reshape(kett_mvmr_harm, timevar="id.exposure", idvar=c("SNP"), direction="wide")
test=as.matrix(test)
col=as.numeric(grep("exposure|SNP", colnames(test)))
test=test[,col]
length(unique(test[,1]))
snps=unique(kett_mvmr_harm$SNP)
effect_alleles=as.data.frame(matrix(1:length(snps),nrow=length(snps),ncol=3))
colnames(effect_alleles)=c("SNP", "eaf.outcome", "remove")
for (i in 1:length(snps))
{
        effect_alleles[i,]=kett_mvmr_harm[i,c("SNP", "eaf.outcome", "remove")]
}
merged=merge(effect_alleles, test, by="SNP")
outcome_dat=read_outcome_data(snps=merged$SNP,filename="UKBB_birthweight",snp_col="RSID", beta_col = "beta",
                              effect_allele_col = "ea", other_allele_col = "nea" , eaf_col = "eaf", pval_col = "p", phenotype_col = "Birthweight")

beta_x_pyr=as.numeric(merged$`beta.exposure.met-d-Pyruvate`)
se_x_pyr=as.numeric(merged$`se.exposure.met-d-Pyruvate`)
beta_y_pyr=as.numeric(outcome_dat$beta.outcome)
se_y_pyr=as.numeric(outcome_dat$se.outcome)

beta_x_glu=as.numeric(merged$`beta.exposure.met-d-Glucose`)
se_x_glu=as.numeric(merged$`se.exposure.met-d-Glucose`)
beta_y_glu=as.numeric(outcome_dat$beta.outcome)
se_y_glu=as.numeric(outcome_dat$se.outcome)

snps=outcome_dat$SNP

png("snp-pyruvate+glucose", width = 800, height = 800)
plot(beta_x_pyr, beta_y_pyr, pch=1, ylab="SNP-BWT assoc.", xlab = "SNP-metabolite assoc.", cex.lab=1.2,
     cex.axis=1.5, cex=1.5, xlim=c(-0.3, 0.2), ylim=c(-0.15, 0.25))
arrows(beta_x_pyr, beta_y_pyr-(se_y_pyr*1.96), beta_x_pyr,
       beta_y_pyr+(se_y_pyr*1.96), length=0, angle=90, code=3, lty=2)
arrows(beta_x_pyr-(se_x_pyr*1.96), beta_y_pyr, beta_x_pyr+
               (se_x_pyr*1.96), beta_y_pyr, length=0, angle=90, code=3, lty=2)
points(beta_x_glu, beta_y_glu, col="red",pch=2)
arrows(beta_x_glu, beta_y_glu-(se_y_glu*1.96), beta_x_glu,
       beta_y_glu+(se_y_glu*1.96), length=0, angle=90, code=3, lty=2,col="red")
arrows(beta_x_glu-(se_x_glu*1.96), beta_y_glu, beta_x_glu+
               (se_x_glu*1.96), beta_y_glu, length=0, angle=90, code=3, lty=2,col="red")
legend("topleft", 
       legend = c("Pyruvate", "Glucose"), 
       col = c("black","red"), 
       pch = c(1,2), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.1, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))
dev.off()

