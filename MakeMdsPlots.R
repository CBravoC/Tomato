#install.packages("git2r")
library(git2r)
r = revparse_single(input_fname,"HEAD")
hash = sha(r)

#Load library with useful data structures and functions:
library(edgeR)

#Load custom functions
source("StyleFunctions.R")

#Read "raw counts" data:
counts = read.delim(input_fname,quote="",header=T,
                    as.is=T,sep="\t",row.names = "gene_name")

#Create edgeR data structure
description_column = grep("description",colnames(counts))
counts = counts[,-description_column]
big_DGEList=DGEList(counts=counts,remove.zeros = TRUE)
#Observe that rows with zero counts in all columns will be removed

#Configure color scheme for plots
all_sample_colors = getColorsForSampleCodes(rownames(big_DGEList$samples))

#Bar plot showing sequencing coverage by sample
main="Sequencing Depth per Sample Library"
ylab="Counts (millions)"
libsizes=big_DGEList$samples$lib.size/10**6
mindepth=round(min(libsizes))
maxdepth=round(max(libsizes))
names(libsizes)=rownames(big_DGEList$samples)
par(mar=c(8,4.1,4.1,2.1))
barplot(libsizes,col=all_sample_colors,las=2,
        main=main,ylab=ylab,cex.names = 0.8,
        ylim=c(0,maxdepth+10))

#MDS to examine sample relatedness
plotMDS(big_DGEList,col=all_sample_colors,main="All samples")
indexes = c(grep("F",names(counts)),
            grep("V",names(counts)))
big_DGEList=DGEList(counts=counts[,indexes],remove.zeros = T)
sample_colors = getColorsForSampleCodes(names(counts[,indexes]))
plotMDS(big_DGEList,col=all_sample_colors,main="F3H-OX and VF36")

#Define a function to cluster samples by genotype
clusterByGenotype = function(genotype,counts) {
  indexes=sapply(names(counts),
                 function(x){getGenotype(x)==genotype})
  sample_colors = getColorsForSampleCodes(names(counts[,indexes]))
  little_DGEList=DGEList(counts[,indexes],
                         remove.zeros = TRUE)
  display_genotype = ""
  if (genotype == "A") {
    display_genotype = "ARE"
  }
  if (genotype == "F") {
    display_genotype = "F3H-OX"
  }
  if (genotype == "V") {
    display_genotype = "VF36"
  }
  main = paste(display_genotype)
  plotMDS(little_DGEList,col=sample_colors,main=main)
}

#ARE
clusterByGenotype("A",counts)

#F3H-OX
clusterByGenotype("F",counts)


#VF36
clusterByGenotype("V",counts)

#Function to cluster samples from the same time point and genotype
clusterByTimeGenotype = function(timepoint,genotype,counts) {
  indexes=sapply(names(counts),
                 function(x){getGenotype(x)== genotype & getTimePoint(x)==timepoint})
  sample_colors = getColorsForSampleCodes(names(counts[,indexes]))
  little_DGEList=DGEList(counts[,indexes],
                         remove.zeros = TRUE)
  display_genotype = ""
  if (genotype == "A") {
    display_genotype = "ARE"
  }
  if (genotype == "F") {
    display_genotype = "F3H-OX"
  }
  if (genotype == "V") {
    display_genotype = "VF36"
  }
  main = paste(display_genotype,
               "at",
               timepoint,
               "minutes")
  plotMDS(little_DGEList,col=sample_colors,main=main)
}

#Cluster genotypes and time points separately

#Cluster ARE samples
timepoints = unique(sapply(names(counts),getTimePoint))
ncols=2
num_plots = length(timepoints)
nrows = num_plots/ncols
par(mfrow=c(nrows,ncols))
genotype="A"
for (timepoint in timepoints) {
    dge_list = clusterByTimeGenotype(timepoint,genotype,counts)
}

#Cluster F3H samples
par(mfrow=c(nrows,ncols))
genotype="F"
for (timepoint in timepoints) {
    dge_list = clusterByTimeGenotype(timepoint,genotype,counts)
}

#Cluster VF36 samples
par(mfrow=c(nrows,ncols))
genotype="V"
for (timepoint in timepoints) {
    dge_list = clusterByTimeGenotype(timepoint,genotype,counts)
}
