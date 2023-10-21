#Define variables
output_fname_prefix = "results/CvT-compare-"
output_fname_suffix = ".txt"
assemblies = c("SL4","SL5")
methods = c("edgeR","DESeq2")
#Analysis and Results
#Load custom functions for analyzing and visualizing differential expression
source("Common.R")

#Load required edgeR and DESeq2 generated data
dfs=list()
for (method in methods) {
  for (assembly in assemblies) {
    thing = getCvS(method=method,assembly=assembly,
                   standardize=T)
    thing$comp = paste(thing$group1,thing$group2,
                       sep="-")
    key = paste(method,assembly,sep=".")
    row.names(thing)=paste(thing$comp,thing$gene,sep=".")
    dfs[[key]]=thing
  }
}

#FDR cutoff variable "Q" for deciding differential expression
Q = 0.05

#Report the number of test results obtained per comparison:
for (key in names(dfs)) {
  thing = dfs[[key]]
  print("COMPARISON")
  print(key)
  print(table(thing$comp))
}
#Show how many tests had adjusted p value of `r Q` or smaller.
#DE genes found by method, assembly:
for (key in names(dfs)) {
  thing = dfs[[key]]
  print("COMPARISON")
  print(key)
  v = thing$padj<=Q
  print(table(thing$comp[v]))
}

#Compare nominal p values:
assembly_synonyms = c("S_lycopersicum_Sep_2019",
                      "S_lycopersicum_Jun_2022")
names(assembly_synonyms)=assemblies
for (assembly in assemblies) {
  key1 = paste(methods[1],assembly,sep=".")
  key2 = paste(methods[2],assembly,sep=".")
  thing1 = dfs[[key1]]
  thing2 = dfs[[key2]]
  v = intersect(row.names(thing1),row.names(thing2))
  plot(-log10(thing1[v,"pvalue"]),
       -log10(thing2[v,"pvalue"]),
       xlab=paste("-log10 pvalue",methods[1]),
       ylab=paste("-log10 pvalue",methods[2]),
       main=paste(assembly,assembly_synonyms[assembly],sep=" - "),
       pch="x")
}

#Examine correspondence of fold-changes:
for (assembly in assemblies) {
  key1 = paste(methods[1],assembly,sep=".")
  key2 = paste(methods[2],assembly,sep=".")
  thing1 = dfs[[key1]]
  thing2 = dfs[[key2]]
  v = intersect(row.names(thing1),row.names(thing2))
  plot(thing1[v,"logFC"],
       thing2[v,"logFC"],
       xlab=paste("log2 fold-change",methods[1]),
       ylab=paste("log2 fold-change",methods[2]),
       main=paste(assembly,assembly_synonyms[assembly],sep=" - "),
       pch="x")
}

#Combine results and includes all the DE genes from each method from the same genome assembly
for (assembly in assemblies) {
  key1 = paste(methods[1],assembly,sep=".")
  key2 = paste(methods[2],assembly,sep=".")
  thing1 = dfs[[key1]]
  thing1 = thing1[thing1$padj<=Q,]
  thing2 = dfs[[key2]]
  thing2 = thing2[thing2$padj<=Q,]
  thing3 = merge(thing1,thing2,
                by="row.names",
                all.x=T,
                all.y=T,
                suffixes=c(paste0(".",methods[1]),
                           paste0(".",methods[2])))
  if (assembly=="SL4") {
    to_write = thing3[,c(2,3,4,5,6,7,14,15,16,19)]
    names(to_write)[c(1:3,10)]=c("gene","group1","group2","description")
    all.sl4 = to_write
  }
  if (assembly=="SL5") {
    to_write = thing3[,c(2,3,4,5,6,7,15,16,17,20,21)]
    names(to_write)[c(1:3,10,11)]=c("gene","group1","group2","description","SL4")
    all.sl5 = to_write
  }
  output_fname=paste0(output_fname_prefix,assembly,output_fname_suffix)
  write.table(to_write,file=output_fname,sep="\t",quote=F,
            row.names = F)
}

#Columns included:

#* gene - the gene measured
#* group1 - the group tested, the control
#* group 2 - the treatment group testing, the treatment
#* padj, pvalue, logFC - adjusted p value, nominal pvalue, log2 fold-change for each method, as indicated in the column name (e.g., padj.edgeR or padj.DESeq2)
#* description - gene description from SL4, if available; NA if not
#* SL4 - only in the SL5 spreadsheet; the SL4 gene name if available

#edgeR and DESeq2 produced similar results

# Session Info
sessionInfo()



