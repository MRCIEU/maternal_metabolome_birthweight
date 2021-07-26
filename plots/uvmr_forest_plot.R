## Forestplot
setwd("/Volumes/MRC-IEU-research/projects/ieu2/p6/106/working/data/data_v2/ukbb_snps")
devtools::install_github("NightingaleHealth/ggforestplot")
library(readr); library(ggplot2); library(data.table);library(TwoSampleMR);library(ggforestplot);library(readxl)

nmr_metabolites=read_excel("nmr_metabolites_20210610.xlsx")
nmr_metabolites=nmr_metabolites[which(nmr_metabolites$Include=="yes"),]

# main analysis
ukbb_uvmr_results=read_csv("naive_uvmr_results.csv")[,-1]

ukbb_uvmr_results$group=NA
for (i in nmr_metabolites$`Biomarker name`)
{
  if (i%in%ukbb_uvmr_results$trait)
  {
    ukbb_uvmr_results$group[which(ukbb_uvmr_results$trait==i)]=nmr_metabolites$Group[which(nmr_metabolites$`Biomarker name`==i)]
  }
}
ukbb_uvmr_results$`Effect estimate (95% CI)`=ukbb_uvmr_results$b
ukbb_uvmr_results=ukbb_uvmr_results[order(ukbb_uvmr_results$trait),]

png("forest_plot_main.png", width = 1000, height = 1400)
ggforestplot::forestplot(df = ukbb_uvmr_results, name = trait, estimate = `Effect estimate (95% CI)`, se = se, pvalue = pval)+
  ggforce::facet_col(facets=~group, scales="free_y", space="free")
dev.off()

# sensitivity forest plot

ket_uvmr_results=read_csv("naive_uvmr_results.csv")[,-1]
multi_uvmr_results=read_csv("naive_uvmr_results.csv")[,-1]

ket_uvmr_results$Dataset="Kettunen"
multi_uvmr_results$Dataset="Multi GWAS"

old=c("18:2, linoleic acid", "22:6, docosahexaenoic acid" ,"Omega-7, omega-9 and saturated fatty acids" ,"3-hydroxybutyrate" ,"Apolipoprotein A-I" ,"Concentration of chylomicrons and largest VLDL particles", 
      "Mono-unsaturated fatty acids" ,"Phosphatidylcholine and other cholines", "Serum total triglycerides", "Total phosphoglycerides", "Free cholesterol", "Serum total cholesterol",
      "Total lipids in chylomicrons and largest VLDL particles")
new=c("Linoleic acid", "Docosahexaenoic acid", "Saturated fatty acids", "3-Hydroxybutyrate", "Apolipoprotein A1", "Concentration of chylomicrons and extremely large VLDL particles",
      "Monounsaturated fatty acids", "Phosphatidylcholines", "Total triglycerides", "Phosphoglycerides", "Total free cholesterol", "Total cholesterol", "Total lipids in chylomicrons and extremely large VLDL")

for (i in 1:length(old))
{
  ket_uvmr_results$trait[grep(paste(old[i]),ket_uvmr_results$trait)]=new[i]
}

old_2=c("18:2, linoleic acid", "22:6, docosahexaenoic acid" ,"Omega-7, omega-9 and saturated fatty acids" ,"3-hydroxybutyrate" ,"Apolipoprotein A-I" , 
        "Mono-unsaturated fatty acids" ,"Phosphatidylcholine and other cholines", "Serum total triglycerides", "Total phosphoglycerides", "Free cholesterol", "Serum total cholesterol",
        "Total lipids in chylomicrons and largest VLDL particles", "Fasting blood glucose", "Triglycerides")
new_2=c("Linoleic acid", "Docosahexaenoic acid", "Saturated fatty acids", "3-Hydroxybutyrate", "Apolipoprotein A1",
        "Monounsaturated fatty acids", "Phosphatidylcholines", "Total triglycerides", "Phosphoglycerides", "Total free cholesterol", "Total cholesterol",
        "Total lipids in chylomicrons and extremely large VLDL", "Glucose", "Total triglycerides")

for (i in 1:length(old_2))
{
  multi_uvmr_results$trait[grep(paste(old_2[i]), multi_uvmr_results$trait)]=new_2[i]
}

#ao=available_outcomes()
ao_1=ao[grep("met-d-", ao$id),]
ao_1=ao_1[which(ao_1$trait%in%ukbb_uvmr_results$id),]

ket_uvmr_results=ket_uvmr_results[ , order(names(ket_uvmr_results))]
multi_uvmr_results=multi_uvmr_results[ , order(names(multi_uvmr_results))]

dataframe=rbind(ket_uvmr_results, multi_uvmr_results)
dataframe$group=NA
for (i in nmr_metabolites$`Biomarker name`)
{
  if (i%in%dataframe$trait)
  {
    dataframe$group[which(dataframe$trait==i)]=nmr_metabolites$Group[which(nmr_metabolites$`Biomarker name`==i)]
  }
}

# manually add groups

dataframe$group[which(dataframe$trait%in%c("Glycated hemoglobin levels", "HDL cholesterol","LDL cholesterol", "Serum creatinine"))]=
                c("Glycolysis related metabolites", "Fluid balance", "Lipoprotein subclasses", "Lipoprotein subclasses")

colnames(dataframe)[1]="Effect estimate (95% CI)"
dataframe$trait[which(is.na(dataframe$group))]
dataframe$group[which(is.na(dataframe$group))]=nmr_metabolites$Group[which(nmr_metabolites$`Biomarker name`==
                                                                             dataframe$trait[which(is.na(dataframe$group))])]
dataframe=dataframe[order(dataframe$trait),]
write_csv(dataframe,"uvmr_kett_multi.csv")
dataframe=read_csv("uvmr_kett_multi.csv")
png("forest_plot_SA.png", width = 1000, height = 1400)
ggforestplot::forestplot(df=dataframe, estimate=`Effect estimate (95% CI)`, pvalue = pval, name=trait, 
                         colour=Dataset, shape=Dataset)+ggforce::facet_col(facets=~group, scales="free_y", space="free")
dev.off()
