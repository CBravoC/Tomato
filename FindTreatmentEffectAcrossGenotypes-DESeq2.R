#Define variables
assemblies=c("SL4","SL5")
Q = 0.05
lfcThreshold=1
outfname_prefix="results/MvW-temp"
#Load functions
source("Common.R")

#Read counts data for assemblies
dfs = list()
for (assembly in assemblies) {
  dfs[[assembly]]=getCounts(assembly=assembly,
                            keep_description = T)
}

#Fit model comparing *are* to VF36 at specific timepoints and treatment
getDEgenes_time <- function(minute) {
  ddss = list()
  for (assembly in assemblies) {
    d = dfs[[assembly]]
    desc_index = which(names(d)=="description")
    d = d[,-desc_index]
    toks = strsplit(names(d),"\\.")
    genotype=sapply(toks,function(x){x[[1]]})
    temperature=sapply(toks,function(x){x[[2]]})
    time=sapply(toks,function(x){x[[3]]})
    v = genotype%in%c("A","V")&time==minute
    d = d[,v]
    coldata=data.frame(genotype=factor(genotype[v],levels=c("V","A")),
                       temperature=factor(temperature[v], levels = c("28","34"))) 
    row.names(coldata)=names(d)
    cts = round(as.matrix(d))
    dds = DESeqDataSetFromMatrix(countData=cts,
                                 colData=coldata,
                                 design=~genotype + temperature + genotype:temperature) 
    featureData = data.frame(gene=rownames(cts))
    mcols(dds) = DataFrame(mcols(dds), featureData)
    dds = DESeq(dds, minReplicatesForReplace=Inf)
    ddss[[assembly]]=dds
  }
  return(ddss)
}

#Collect results function
results_time <- function(ddss, minute) {
  rss=list()
  for (assembly in assemblies) {
    dds=ddss[[assembly]]
    rs = results(dds,alpha=Q,lfcThreshold = lfcThreshold,
                 name="genotypeA.temperature34") 
    rs = rs[!is.na(rs$padj),]
    rs=cbind(gene=row.names(rs),rs)
    rs["time"] = minute 
    rs$description = dfs[[assembly]][rs$gene,"description"]
    if (assembly == "SL5") {
      SL4_gene_names = getSL4GeneNames(rs$description)
      rs$SL4 = SL4_gene_names
    }
    o = order(rs$pvalue)
    rs = rs[o,]
    rss[[assembly]]=data.frame(rs)
  }
  return(rss)
}

#15 minute Results
min15 <- getDEgenes_time(15)
res15 <- results_time(min15, 15)

#30 minutes Results
min30 <- getDEgenes_time(30)
res30 <- results_time(min30, 30)

#45 minutes Results
min45 <- getDEgenes_time(45)
res45 <- results_time(min45, 45)

#75 minutes Results 
min75 <- getDEgenes_time(75)
res75 <- results_time(min75, 75)

#Combine the dataframes 

combined_results <- rbind(res15, res30, res45, res75)
SL4 <- combined_results[,1]
SL5 <- combined_results[,2]

#View and Save the Results by assembly, write a file with results
library(openxlsx)
assembly="SL4"
outfname=paste0(outfname_prefix,"-",assembly,".xlsx")
data_frames <- list("15 minutes" = SL4$res15, "30 minutes" = SL4$res30, 
                    "45 minutes" = SL4$res45,"75 minutes" = SL4$res75 )
write.xlsx(data_frames, file = outfname)

assembly="SL5"
outfname=paste0(outfname_prefix,"-",assembly,".xlsx")
data_frames <- list("15 minutes" = SL5$res15, "30 minutes" = SL5$res30, 
                    "45 minutes" = SL5$res45,"75 minutes" = SL5$res75 )
write.xlsx(data_frames, file = outfname)

#Session Info
sessionInfo()