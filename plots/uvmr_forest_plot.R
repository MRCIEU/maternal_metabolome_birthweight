## Forestplot

kett_uvmr_results=read_csv("kett_uvmr_results.csv")
ukbb_uvmr_results=read_csv("ukbb_uvmr_results.csv")
kett_uvmr_results$Dataset="Kettunen"
ukbb_uvmr_results$Dataset="UKBB"
nmr_metabolites_UKBB=read_excel("nmr_metabolites.xlsx")
ao=available_outcomes()
ao_1=ao[grep("met-d-", ao$id),]
ao_1=ao_1[which(ao_1$trait%in%ukbb_uvmr_results$id),]

colnames(ukbb_uvmr_results)[c(1,11)]=c("trait","id")
# fix mismatching names
old=c("18:2, linoleic acid", "22:6, docosahexaenoic acid" ,"Omega-7, omega-9 and saturated fatty acids" ,"3-hydroxybutyrate" ,"Apolipoprotein A-I" ,"Concentration of chylomicrons and largest VLDL particles", 
      "Mono-unsaturated fatty acids" ,"Phosphatidylcholine and other cholines", "Serum total triglycerides", "Total phosphoglycerides", "Free cholesterol", "Serum total cholesterol",
      "Total lipids in chylomicrons and largest VLDL particles")
new=c("Linoleic acid", "Docosahexaenoic acid", "Saturated fatty acids", "3-Hydroxybutyrate", "Apolipoprotein A1", "Concentration of chylomicrons and extremely large VLDL particles",
      "Monounsaturated fatty acids", "Phosphatidylcholines", "Total triglycerides", "Phosphoglycerides", "Total free cholesterol", "Total cholesterol", "Total lipids in chylomicrons and extremely large VLDL")

for (i in 1:length(old))
{
  kett_uvmr_results$trait[grep(paste(old[i]),kett_uvmr_results$trait)]=new[i]
}

ukbb_uvmr_results=ukbb_uvmr_results[,which(colnames(ukbb_uvmr_results)%in%colnames(kett_uvmr_results))]
kett_uvmr_results=kett_uvmr_results[,which(colnames(kett_uvmr_results)%in%colnames(ukbb_uvmr_results))]

ukbb_uvmr_results=ukbb_uvmr_results[ , order(names(ukbb_uvmr_results))]
kett_uvmr_results=kett_uvmr_results[ , order(names(kett_uvmr_results))]

dataframe=rbind(ukbb_uvmr_results, kett_uvmr_results)
dataframe$group=NA
for (i in nmr_metabolites$`Biomarker name`)
{
  if (i%in%dataframe$trait)
  {
    dataframe$group[which(dataframe$trait==i)]=nmr_metabolites$Group[which(nmr_metabolites$`Biomarker name`==i)]
  }
}

names(dataframe)[1]="Effect estimate (95% CI)"
dataframe$trait[which(is.na(dataframe$group))]
dataframe$group[which(is.na(dataframe$group))]=nmr_metabolites$Group[which(nmr_metabolites$`Biomarker name`==
                                                                             dataframe$trait[which(is.na(dataframe$group))])]
dataframe=dataframe[order(dataframe$trait),]

png("forest_plot.png", width = 1000, height = 1400)
ggforestplot::forestplot(df=dataframe, estimate=`Effect estimate (95% CI)`, pvalue = pval, name=trait, 
                         colour=Dataset, shape=Dataset)+ggforce::facet_col(facets=~group, scales="free_y", space="free")
dev.off()
