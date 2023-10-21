#Functions to be called later from DEscripts
#Load libraries
library(stringr)
library(readr)
library(readxl)
library(edgeR)
library(DESeq2)
library(EnhancedVolcano)

#Filter significant data
prettyNumber = function(x,signif_digits=NA,
                        big.mark=",") {
  if (!is.na(signif_digits)) {
    to_return = format(signif(x,signif_digits),
                       big.mark=big.mark)
  }
  else {
    to_return = format(x,big.mark=big.mark)
  }
  return(to_return)
}

getSampleExcelSpreadsheet = function(fname="Documentation/muday-144_sample_sheet.xlsx") {
  sample_sheet = read_excel(fname)
  return(sample_sheet)
}

getColorsForSampleCodes = function(sample_names,
                                   sample_sheet_fname=NULL,
                                   sample_sheet=NULL) {
  if (is.null(sample_sheet)) {
    if (is.null(sample_sheet_fname)) {
      sample_sheet = getSampleExcelSpreadsheet()      
    }
    else {
      sample_sheet = getSampleExcelSpreadsheet(sample_sheet_fname)
    }
  }
  sample_codes = sample_sheet$`Sample Code` 
  sample_colors = sample_sheet$Color
  names(sample_colors) = sample_codes
  to_return = sample_colors[sample_names]
  return(to_return)
}

#Need to test this
getSampleColors = function(sample_sheet_fname=NULL,
                           sample_sheet=NULL,
                           hex=F) {
  if (is.null(sample_sheet)) {
    if (!is.null(sample_sheet_fname)) {
      sample_sheet = getSampleExcelSpreadsheet(fname=sample_sheet_fname)
    } 
    else {
      sample_sheet=getSampleExcelSpreadsheet()
    }
  }
  sample_colors = sample_sheet$Color
  sample_codes = sample_sheet$`Sample Code`
  names(sample_colors)=sample_codes
  if (hex) {
#https://gist.github.com/mbannert/e9fcfa86de3b06068c83
    rgbs=col2rgb(sample_colors)
    hexes=apply(rgbs,2,
                function(v){rgb(v[1],v[2],v[3], maxColorValue = 255)})
    return(hexes)
  }
  else {
    return(sample_colors)
  }
}

getSampleGroups = function(sample_names) {
  sample_groups=sub("\\.[789]$",sample_names,replacement="")
  return(sample_groups)
}

#Example input: "F.28.15.7"
getGenotype = function(sample_name) {
  value=strsplit2(sample_name,"\\.")[1,1]
  return(value)
}

getTreatment=function(sample_name) {
  value=strsplit2(sample_name,"\\.")[1,2]
  ifelse(value=='28',"control","heat stress")
}

getTimePoint=function(sample_name) {
  value = getTreatmentDuration(sample_name)
  return(value)
}

getTreatmentDuration=function(sample_name) {
  value = strsplit2(sample_name,"\\.")[1,3]
  return(value)
}

getReplicate=function(sample_name) {
  return (strsplit2(sample_name,"\\.")[1,4])
}

getVarieties = function(sample_sheet_fname="Documentation/muday-144_sample_sheet.xlsx") {
  sample_sheet=getSampleExcelSpreadsheet(sample_sheet_fname)
  sample_codes = sample_sheet$`Sample Code`
  x = sapply(sample_names,getVariety)
  varieties = unique(x)
  return(varieties)
}

getVarietyCodes = function() {
  varieties = getVarieties()
  varietyCodes = sapply(varieties,function(x){substr(x,1,1)})
  return(varietyCodes)
}

#Need to update this
getTreatmentCodes = function() {
  return(c("C","S"))
}

#Get counts per gene per sample
#Return data frame of counts with gene names as row names
#Keep description column if keep_description is TRUE
#Default is don't keep it
getCounts = function(counts_fname=NULL,
                     keep_description=FALSE,
                     assembly="SL5") {
  if (is.null(counts_fname)) {
    counts_fname=paste0("results/muday-144-",assembly,
                        "_counts-salmon.txt")
  }  
  data=read.delim(counts_fname,
                  comment.char = "#",
                  row.names = "gene_name") 
  if (!keep_description) {
    data = data[,1:ncol(data)-1] 
  }
  return(data)
}

#Get scaled counts per gene per sample
#Return data frame of scaled counts with gene names as row names
#Keep description column if keep_description is TRUE
#Default is don't
getScaledCounts = function(counts_fname="results/muday-144-SL5_counts-salmon_scaled.txt",
                           keep_description=FALSE) {
  data=read.delim(counts_fname,
                  comment.char = "#",
                  row.names = "gene_name") 
  if (!keep_description) {
    data = data[,1:ncol(data)-1] 
  }
  return(data)
}

#Get treatment versus control differential expression 
#Result for simple C versus S comparison
#If standardize, change some column names as shown below
getCvS = function(method="edgeR",assembly="SL5",standardize=F) {
  if (!method %in% c("edgeR","DESeq2") | !assembly %in% c("SL5","SL4")) {
    return(NULL)
  }
  else {
    fname = paste("CvT",method,assembly,sep="-")
    fname = paste(fname,"txt",sep=".")
    fname = paste("results",fname,sep="/")
    data=read.delim(fname,comment.char = "#",header=T,quote="",as.is=T)
    if (standardize) {
      if (method=="edgeR") {
        index = which(names(data)=="PValue")
        names(data)[index]="pvalue"
        index = which(names(data)=="Q")
        names(data)[index]="padj"
      }
      if (method=="DESeq2") {
        index = which(names(data)=="log2FoldChange")
        names(data)[index]="logFC"
      }
    }
    return(data)
  }
}

#Need to test this
#Retrieve a named character vector where names are
#Gene names and values are IGB links
getIgbLinks=function(assembly="SL5") {
  v=c("SL5","S_lycopersicum_Jun_2022",
      "SL4","S_lycopersicum_Sep_2019")
  lookup_assembly=c(rep(v[2],2),rep(v[4],2))
  names(lookup_assembly)=v
  annots_fname=paste(lookup_assembly[assembly],
                     "bed.gz",sep=".")
  annots_path=paste("..","..","quickload",
                    annots_fname,sep="/")
  annots=read.delim(annots_path,header=F,as.is=T,
                    quote="",sep="\t")
  annots=annots[,c(1,2,3,4,6,13)]
  names(annots)=c('seq','start','end','tx_id','strand',
                  'gene_id')
  o=order(annots$seq,annots$start)
  annots=annots[o,]
  minstarts=annots[!duplicated(annots$gene_id),
                   c('gene_id','seq','start')]
  o=order(annots$seq,annots$end,decreasing=T)
  annots=annots[o,]
  maxends=annots[!duplicated(annots$gene_id),
                 c('gene_id','end','strand')]
  generegions=merge(minstarts,maxends,by='gene_id')
  links=paste0('http://localhost:7085/UnibrowControl?version=&seqid=',
               generegions$seq,'&start=',generegions$start-200,'&end=',generegions$end+200)
  names(links)=generegions$gene_id
  return(links)
}

#Need to update this
#Get a data frame containing scaled counts as 
#Fragments per million
#fname - data file with counts
#ave - if T, include average per sample group
#sd - if T, include standard deviation per sample group
#returns data frame where rownames are gene ids, e.g., "Solyc00g005280.1"
#Gene name and gene description are included as columns (not factors)
#If averages or sd are requested, output includes average or sd for each sample type
getFPM=function(fname="../Count/results/fpm.txt",
                ave=F,sd=F){
  row.names(d)=d$gene
  if (!ave) {
    d=d[,-grep("ave",names(d))]
  }
  if (!sd) {
    d=d[,-grep("sd",names(d))]
  }
  return(d)
}

#Need to test this
getSampleNames = function(counts_fname) {
  d = getCounts(counts_fname)
  sample_names = colnames(d)
  return(sample_names)
}

#Need to test this
getVariety = function(sample_name) {
  letter_code=substr(sample_name,1,1)
  if (letter_code=='V') {
    return("VF36")
  }
  if (letter_code=="F") {
    return("F3H-OX3")
  }
  if (letter_code=="A") {
    return("are")
  }
  else {
    return(sample_name)
  }
}

#Need to test this
getTreatment=function(sample_name) {
  treatment=substr(sample_name,3,4)
  return(treatment)
}

#Need to test this
getReplicate=function(sample_name) {
  letter_code=substr(sample_name,9,10)
  return(letter_code)
}

#Need to test this
getDescription=function(sample_name) {
  descr=paste(sample_name,"-",getVariety(sample_name),getTreatment(sample_name),
              " degC",
              "replicate",getReplicate(sample_name))
  return(descr)
}

#Need to test this
getCodeBook = function(counts_fname=NULL,
                       sample_sheet_fname=NULL) {
  
  samples = names(getCounts(counts_fname))
  group=str_remove_all(samples,".\\d")
  filename=rep("a-file-name.fastq.gz",length(samples))
  variety = sapply(samples,getVariety)
  treatment=sapply(samples,getTreatment)
  replicate=sapply(samples,getReplicate)
  col=getSampleColors(sample_sheet_fname=sample_sheet_fname,
                      hex=T)
  col=col[samples]
  description=sapply(samples,getDescription)
  codebook = data.frame(sample=samples,group,filename,
                        variety,treatment,replicate,
                        color=col,
                        description)
  return(codebook)
}

getGeneDescription = function(counts_fname) {
  d = read.delim(counts_fname,quote="",header=F,
                 as.is=T,sep="\t")
  d = d[,c(1,ncol(d))]
  names(d)=c("gene","description")
  row.names(d)=d$gene
  return(d)
}


#Make a barplot showing expression in each sample for the
#Given gene dat requires: columns for each sample in the experiment description column
makeBarPlot=function(gene,dat,ylab,main=NA,beside=F) {
  description_index = which(names(dat)=="description")
  if (length(description_index)==1) {
    description_index = description_index[1]
    if (is.na(main)) {
      main = dat[gene,description_index]
    }
    dat = dat[,-(description_index)]
  }
  sample_colors = getSampleColors()[names(dat)]
  rowdata=dat[gene,names(dat)]
  rowdata=as.numeric(rowdata)
  names(rowdata)=names(sample_colors)
  if (is.na(main)) {
    main = gene
  }
  barplot(rowdata,
          beside=beside,
          main=main,
          col=sample_colors,
          las=2,
          ylab=ylab)
}

getSL4GeneNames = function(sl5_description) {
  names_matrix = strsplit2(sl5_description,"ITAG4.0:")
  to_return = names_matrix[,2]
  boolean_vec = to_return == ""
  to_return[boolean_vec]=names_matrix[boolean_vec,1]
  return(to_return)
}

getDeGenes = function(counts,group1_name,group2_name,
                      method="DESeq2") {
  to_return = NULL
  if (method=="DESeq2") {
    to_return = compareTwoGroupsDESeq2(counts,group1_name,
                                       group2_name) 
  }
  if (method=="edgeR") {
    to_return = compareTwoGroupsEdgeR(counts,group1_name,
                                      group2_name)
  }
  return(to_return)
}

compareTwoGroupsDESeq2 = function(counts,group1_name,group2_name) {
  group1_indexes=grep(group1_name,names(counts)) 
  group2_indexes=grep(group2_name,names(counts))
  indexes=c(group1_indexes,group2_indexes)
  
  condition = factor(c(rep(group1_name,length(group1_indexes)),
                       rep(group2_name,length(group2_indexes))))
  m = round(as.matrix(counts[,indexes]))
  dds = DESeqDataSetFromMatrix(m,
                               DataFrame(condition),
                               ~ condition)
  featureData = data.frame(gene=rownames(m))
  mcols(dds) = DataFrame(mcols(dds),featureData)
  dds = DESeq(dds, minReplicatesForReplace=Inf)
  results_table = results(dds,
                          alpha=0.05)
  results_table = results_table[!is.na(results_table$padj),]
  results_table$gene=row.names(results_table)
  results_table$group1=group1_name
  results_table$group2=group2_name
  results_table$description = counts[results_table$gene,
                                     "description"]
  o = order(results_table$pvalue)
  results_table = results_table[o,]
  d = data.frame(results_table)
  return(d)
}


compareTwoGroupsEdgeR  = function(counts,group1_name,group2_name) {
  group1_indexes=grep(group1_name,names(counts)) 
  group2_indexes=grep(group2_name,names(counts))
  indexes=c(group1_indexes,group2_indexes)
  group = factor(c(rep(group1_name,length(group1_indexes)),
                   rep(group2_name,length(group2_indexes))))
  little_DGEList = DGEList(counts=counts[,indexes],group=group)
  keep = filterByExpr(little_DGEList)
  little_DGEList = little_DGEList[keep,,keep.lib.sizes=FALSE]
  little_DGEList = calcNormFactors(little_DGEList)
  design=model.matrix(~group)
  little_DGEList = estimateDisp(little_DGEList,design)
  fit = glmFit(little_DGEList, design)
  lrt=glmLRT(fit,coef=2)
  results=lrt$table
  results$Q=p.adjust(results$PValue,method="fdr")
  results$gene=row.names(results)
  results$group1=group1_name
  results$group2=group2_name
  results$description = counts[results$gene,"description"]
  o = order(results$Q)
  results = results[o,]
}

volcano_plot<- function(d, title, method="DESeq2",Q=0.05) {
  colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  if (method=="DESeq2") {
    lfc_colname = "log2FoldChange"
    padj_colname = "padj"
    p_colname = "pvalue"
  }
  if (method=="edgeR") {
    lfc_colname = "logFC"
    padj_colname = "Q"
    p_colname = "PValue"
  }
  pCutoffCol = padj_colname
  if (!any(d[,pCutoffCol]<=Q)) {
    Q = round(min(d$padj),2)+0.01
  }
  lfcmin = min(d[,lfc_colname])
  lfcmax = max(d[,lfc_colname])
  #x_max = -log10(min(d[,p_colname]))
  keyvals <- ifelse(
    d[,lfc_colname] < -1, '#56B4E9',
    ifelse(d[,lfc_colname] > 1, '#CC79A7','black'))
  names(keyvals)[keyvals == '#CC79A7'] = 'Up-regulated'
  names(keyvals)[keyvals == 'black'] = 'Mid-regulated'
  names(keyvals)[keyvals == '#56B4E9'] = 'Down-regulated'
  v = EnhancedVolcano(d,
                      lab = d$gene,
                      x = lfc_colname,
                      y = p_colname,
                      selectLab = d[which(names(keyvals) %in% c('Up-regulated', 'Down-regulated')),
                                    "gene"],
                      colCustom = keyvals,
                      colAlpha = 1,
                      shape = c(4, 1, 6, 3),
                      pCutoff = Q,  
                      FCcutoff = 1.0,
                      pCutoffCol = pCutoffCol, # 
                      gridlines.minor=FALSE, 
                      gridlines.major=FALSE,
                      col = colorBlindGrey8,
                      #xlim = c(-10, 10),
                      #ylim = c(0,pmax),
                      #ylab = "-log(p)",
                      title = title, 
                      subtitle = paste("Horizontal dashed line shows FDR",Q),
                      legendPosition = "right")
  return(v)
}


