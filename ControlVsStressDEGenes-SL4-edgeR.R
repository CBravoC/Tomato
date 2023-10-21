#Define variables
genome_version="SL4"
counts_fname = paste0("results/muday-144-",genome_version,"_counts-salmon.txt")
sample_sheet_fname = "Documentation/muday-144_sample_sheet.xlsx"
all = NULL
library(git2r)
r = revparse_single(counts_fname,"HEAD")
hash = sha(r)
method="edgeR"
out_fname = paste0("results/CvT-",method,"-",genome_version,".txt")

#Load custom functions:
source("Common.R")

#Load required counts data:
counts=getCounts(counts_fname,keep_description = T)

# Compute and tabulate differential expression
Q=0.05

#*are* at 15min
group1_name = "A.28.15"
group2_name = "A.34.15"
results_table = getDeGenes(counts,group1_name,group2_name,method=method)
volcano_A1 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = results_table

#*are* at 30min
group1_name = "A.28.30"
group2_name = "A.34.30"
results_table = getDeGenes(counts,group1_name,group2_name,method=method)
volcano_A2 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#*are* at 45min
group1_name = "A.28.45"
group2_name = "A.34.45"
results_table = getDeGenes(counts,group1_name,group2_name,method=method)
volcano_A3 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#*are* at 75min
group1_name = "A.28.75"
group2_name = "A.34.75"
results_table = getDeGenes(counts,group1_name,group2_name,method=method)
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

#VF36 at 15min

group1_name = "V.28.15"
group2_name = "V.34.15"
results_table = getDeGenes(counts,group1_name,group2_name,method=method)
volcano_A1 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#VF36 at 30min
group1_name = "V.28.30"
group2_name = "V.34.30"
results_table = getDeGenes(counts,group1_name,group2_name,method=method)
volcano_A2 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#VF36 at 45min
group1_name = "V.28.45"
group2_name = "V.34.45"
results_table = getDeGenes(counts,group1_name,group2_name,method=method)
volcano_A3 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#VF36 at 75min
group1_name = "V.28.75"
group2_name = "V.34.75"
results_table = getDeGenes(counts,group1_name,group2_name,method=method)
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
results_table = getDeGenes(counts,group1_name,group2_name,method=method)
volcano_A1 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#F3H-OX3, F3H overexpression at 30min
group1_name = "F.28.30"
group2_name = "F.34.30"
results_table = getDeGenes(counts,group1_name,group2_name,method=method)
volcano_A2 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#F3H-OX3, F3H overexpression at 45min
group1_name = "F.28.45"
group2_name = "F.34.45"
results_table = getDeGenes(counts,group1_name,group2_name,method=method)
volcano_A3 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#F3H-OX3, F3H overexpression at 75min
group1_name = "F.28.75"
group2_name = "F.34.75"
results_table = getDeGenes(counts,group1_name,group2_name,method=method)
volcano_A4 <- volcano_plot(results_table, paste(group1_name,"vs",group2_name),
                           method=method)
all = rbind(all,results_table)

#Volcano plots summarizing the above results_table
A_combined<- ggpubr::ggarrange(volcano_A1, volcano_A2, volcano_A3, volcano_A4, # list of plots
                  labels = "AUTO",
                  font.label = list(size = 30), 
                  common.legend = T, # COMMON LEGEND
                  legend = "right", # legend position
                  align = "hv", # Align them both, horizontal and vertical
                  nrow = 2, 
                  ncol = 2)
A_combined

#Format and organize results:
all$logFC=round(all$logFC,3)
all$logCPM=round(all$logCPM,3)
all$Q=signif(all$Q,3)
all$PValue=signif(all$PValue,3)
all = all[,c("gene","group1","group2","Q","PValue","logFC",
             "logCPM","description")]

#Write all the results to a data file
write.table(all,file=out_fname,quote=F,row.names = F,sep="\t")

#Explanation of columns:
#* gene - `r genome_version` gene measured
#* group 1 - control group
#* group 2 - treatment group
#* Q - false discovery rate; adjusted p-value computed using method of Benjamini and Hochberg
#* PValue - nominal (unadjusted) p value 
#* logFC - log2(group 2 average/group 1 average)
#* logCPM - ?
#* description - gene description 

#Session info
sessionInfo()

