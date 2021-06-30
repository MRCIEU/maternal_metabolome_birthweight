## Heatmaps
library(readr);library(gplots)

# I. UKBB
mvmr_harm=read_csv("ukbb_mvmr_harm")
setDT(mvmr_harm)
test=reshape(mvmr_harm, timevar="id.exposure", idvar=c("SNP"), direction="wide")
test=as.matrix(test)
test=as.data.frame(test)[,grep("beta.exposure", colnames(test))]
ao=available_outcomes()
ao_1=ao[which(ao$year==2020&ao$author=="Borges CM"),]
colnames(test)=sub("beta.exposure.","", colnames(test))

test=apply(test, 2, function(x) as.numeric(as.character(x)))
sapply(test, class)
for (j in 1:nrow(ao_1))
{
  colnames(test)=gsub(ao_1$id[j], ao_1$trait[j], colnames(test))
  colnames(test)=gsub(ao$id[j], ao_1$trait[j], colnames(test))
}
cormatrix=cor(test)
pdf('heatmap_ukbb_snp', 12, 12)
heatmap.2(cormatrix, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram='none',
          Rowv=TRUE, Colv=TRUE, keysize=1, margins=c(18.5,20))
graphics.off()

# II. Kett
mvmr_harm=read_csv("kett_snps_mvmr_harm")
setDT(mvmr_harm)
test=reshape(mvmr_harm, timevar="id.exposure", idvar=c("SNP"), direction="wide")
test=as.data.frame(test)[,grep("beta.exposure", colnames(test))]
ao=available_outcomes()
ao_1=ao[which(ao$year==2020&ao$author=="Borges CM"),]
colnames(test)=sub("beta.exposure.","", colnames(test))

test=apply(test, 2, function(x) as.numeric(as.character(x)))
sapply(test, class)
for (j in 1:nrow(ao_1))
{
  colnames(test)=gsub(ao_1$id[j], ao_1$trait[j], colnames(test))
  colnames(test)=gsub(ao$id[j], ao_1$trait[j], colnames(test))
}
cormatrix=cor(test)
pdf('heatmap_kett_snp', 12, 12)
heatmap.2(cormatrix, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram='none',
          Rowv=TRUE, Colv=TRUE, keysize=1, margins=c(18.5,20))
graphics.off()

# III. Multi
mvmr_harm=read_csv("multi_snps_mvmr_harm")
setDT(mvmr_harm)
test=reshape(mvmr_harm, timevar="id.exposure", idvar=c("SNP"), direction="wide")
test=as.matrix(test)
test=as.data.frame(test)[,grep("beta.exposure", colnames(test))]
ao=available_outcomes()
ao_1=ao[which(ao$year==2020&ao$author=="Borges CM"),]
colnames(test)=sub("beta.exposure.","", colnames(test))

test=apply(test, 2, function(x) as.numeric(as.character(x)))
sapply(test, class)
for (j in 1:nrow(ao_1))
{
  colnames(test)=gsub(ao_1$id[j], ao_1$trait[j], colnames(test))
  colnames(test)=gsub(ao$id[j], ao_1$trait[j], colnames(test))
}
cormatrix=cor(test)
pdf('heatmap_multi_snp', 12, 12)
heatmap.2(cormatrix, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram='none',
          Rowv=TRUE, Colv=TRUE, keysize=1, margins=c(18.5,20))
graphics.off()
