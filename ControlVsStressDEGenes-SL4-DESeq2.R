#Load custom functions for analyzing and visualizing differential expression
source("Common.R")

#Load required counts data:
counts=getCounts(counts_fname,keep_description = T)

#Define FDR threshold for deciding whether a gene is DE
Q=0.05

#Genotype *are* at 15min
group1_name = "A.28.15"
group2_name = "A.34.15"
results_table = getDeGenes(counts,group1_name,group2_name)
volcano_A1 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = results_table

#Genotype *are* at 30min
group1_name = "A.28.30"
group2_name = "A.34.30"
results_table = getDeGenes(counts,group1_name,group2_name)
volcano_A2 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#Genotype *are* at 45min
group1_name = "A.28.45"
group2_name = "A.34.45"
results_table = getDeGenes(counts,group1_name,group2_name)
volcano_A3 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#Genotype *are* at 75min
group1_name = "A.28.75"
group2_name = "A.34.75"
results_table = getDeGenes(counts,group1_name,group2_name)
volcano_A4 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#Volcano plots summarizing above results
A_combined<- ggpubr::ggarrange(volcano_A1, volcano_A2, volcano_A3, volcano_A4, #list of plots
                  labels = "AUTO",
                  font.label = list(size = 30), 
                  common.legend = T, # COMMON LEGEND
                  legend = "right", # legend position
                  align = "hv", # Align them both, horizontal and vertical
                  nrow = 2, 
                  ncol = 2)
A_combined

#VF36 at 15min
group1_name = "V.28.15"
group2_name = "V.34.15"
results_table = getDeGenes(counts,group1_name,group2_name)
volcano_A1 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#VF36 AT 30min
group1_name = "V.28.30"
group2_name = "V.34.30"
results_table = getDeGenes(counts,group1_name,group2_name)
volcano_A2 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#VF36 AT 45min
group1_name = "V.28.45"
group2_name = "V.34.45"
results_table = getDeGenes(counts,group1_name,group2_name)
volcano_A3 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#VF36 AT 75min
group1_name = "V.28.75"
group2_name = "V.34.75"
results_table = getDeGenes(counts,group1_name,group2_name)
volcano_A4 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#Volcano plots summarizing the above results
A_combined<- ggpubr::ggarrange(volcano_A1, volcano_A2, volcano_A3, volcano_A4, # list of plots
                  labels = "AUTO",
                  font.label = list(size = 30), 
                  common.legend = T, # COMMON LEGEND
                  legend = "right", # legend position
                  align = "hv", # Align them both, horizontal and vertical
                  nrow = 2, 
                  ncol = 2)
A_combined

#F3H-OX3, F3H overexpression at 15min
group1_name = "F.28.15"
group2_name = "F.34.15"
results_table = getDeGenes(counts,group1_name,group2_name)
volcano_A1 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#F3H-OX3, F3H overexpression at 30min
group1_name = "F.28.30"
group2_name = "F.34.30"
results_table = getDeGenes(counts,group1_name,group2_name)
volcano_A2 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#F3H-OX3, F3H overexpression at 45min
group1_name = "F.28.45"
group2_name = "F.34.45"
results_table = getDeGenes(counts,group1_name,group2_name)
volcano_A3 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#F3H-OX3, F3H overexpression at 75min
group1_name = "F.28.75"
group2_name = "F.34.75"
results_table = getDeGenes(counts,group1_name,group2_name)
volcano_A4 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#Volcano plots summarizing the above results
A_combined<- ggpubr::ggarrange(volcano_A1, volcano_A2, volcano_A3, volcano_A4, # list of plots
                  labels = "AUTO",
                  font.label = list(size = 30), 
                  common.legend = T, # COMMON LEGEND
                  legend = "right", # legend position
                  align = "hv", # Align them both, horizontal and vertical
                  nrow = 2, 
                  ncol = 2)
A_combined

#Organize results and round numeric results to 3 significant digits:
all = all[,c("gene","group1","group2","baseMean","padj","pvalue","log2FoldChange","lfcSE","stat","description")]
for (i in 4:9) {
  all[,i]=signif(all[,i],3)
}

#Write the data file:
write.table(all,file=out_fname,quote=F,row.names = F,sep="\t")

#Explanation of columns:
  
#* gene - r genome_version gene measured
#* group 1 - control group
#* group 2 - treatment group
#* baseMean - mean across samples
#* padj - false discovery rate; adjusted p-value computed using method of Benjamini and Hochberg
#* log2FoldChange - log2(group 2 average/group 1 average)
#* lfcsE - log2FoldChange standard error
#* stat - test statistic used to assess significance 
#* description - gene description 

#Session info
sessionInfo()


