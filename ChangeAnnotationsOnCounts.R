#Define variables
genome_version = "SL4"
annotations_fname = ifelse(genome_version=="SL5",
                           "../ExternalDataSets/S_lycopersicum_Jun_2022.bed.gz",
                           "../ExternalDataSets/S_lycopersicum_Sep_2019.bed.gz")
#Load a library
library(readxl)

#Define variables
data_prefix="muday-144-"
data_suffix="_salmon.merged.gene_counts.tsv"
counts_fname=paste0(data_prefix,genome_version,data_suffix)
sample_sheet = './Documentation/muday-144_sample_sheet.xlsx'
sample_renaming_table = "./Documentation/Muday-lab-RNA-samples-for-sample-name-conversion.xlsx"
output_fname = paste0("./results/",data_prefix,genome_version,"_counts-salmon.txt")

#Read gene counts file
counts=read.delim(counts_fname,header = T,row.names = "gene_id")

#Read gene model annotation file
annotations = read.delim(annotations_fname,quote="",
                         header=F,as.is=T,sep="\t")[,13:14]
num_genes = length(unique(annotations[,1]))

#Merge annotations and counts data
names(annotations)=c("gene_name","description")
new_counts=merge(counts,annotations,
                 by.x="gene_name",
                 by.y="gene_name")

#Create a vector with desired column names
genotype = rep(c(rep("F",8),rep("V",8),rep("A",8)),3)
temperature = rep(c("28","34"),36)
time = rep(c(rep("15",2),rep("30",2),rep("45",2),rep("75",2)),9)
replicate = c(rep("7",24),rep("8",24),rep("9",24))
new_colnames = paste(genotype,temperature,
                     time,replicate,sep=".")
new_colnames = c("gene_name",new_colnames,"description")

#Sample codes and new column names will be
new_colnames

#Definite a function, returns `TRUE` if a proposed column name correctly renames a current column name (renaming due to Azenta's mislabeling of samples)
areSameSample = function(a_row) {
  new = a_row[1]
  old = a_row[2]
  if (new==old) {
    return(TRUE)
  }
  new_vals = unlist(strsplit(new,"\\."))
  new_genotype = new_vals[1]
  new_temperature = new_vals[2]
  new_time = new_vals[3]
  new_replicate = new_vals[4]
  old_vals = unlist(strsplit(old,"\\."))
  old_genotype = old_vals[2]
  old_replicate = old_vals[1]
  old_time = old_vals[3]
  old_temperature = old_vals[5]
  if (new_time != old_time) {
    print(new_time)
    print(old_time)
    return(FALSE)
  }
  if (new_temperature != substr(old_temperature,1,2)) {
    print(new_temperature)
    print(old_temperature)
    return(FALSE)
  }
  if (new_replicate!=substr(old_replicate,2,3)) {
    print(new_replicate)
    print(old_replicate)
    return(FALSE)
  }
  new_to_old_lookup = c("are","OE3","VF36")
  names(new_to_old_lookup)=c("A","F","V")
  if (new_to_old_lookup[new_genotype] !=old_genotype) {
    print(new_to_old_lookup[new_genotype])
    print(old_genotype)
    return(FALSE)
  }
  return(TRUE)
}

#Apply new funciton
current_colnames = colnames(new_counts)
to_test = data.frame(new_colnames,current_colnames)
result = apply(to_test,1,areSameSample)
if (sum(result)==length(result)) {
  colnames(new_counts) = new_colnames
} else {
  print("WARNING: Can't rename column names.")
}

#Sample codes: * [genotype].[temperature].[minutes].[replicate]
genotypes=c("V","F","A")
names(genotypes)=c("VF36","OE3","are")
translate = function(x) {
  x1 = sub("#",'',x)
  x2 = sub("C",'',x1)
  v = strsplit(x2," ")[[1]]
  v1 = c(genotypes[v[2]],v[5],v[3],v[1])
  sample_code=paste(v1,collapse = ".")
  return(sample_code)
}
key = read_excel(sample_renaming_table,skip=3)
key=subset(key,!is.na(key[,1]) & !is.na(key[,3]))

#Use Muday lab sample renaming key to rename `new_counts` columns:
original=sapply(unlist(as.vector(key[,2])),translate)
desired = sapply(unlist(as.vector(key[,3])),translate)
names(desired)=original
current_names = names(new_counts)
corrected_names = rep(NA,length(current_names))
for (iter in seq(1,length(current_names))) {
  name = current_names[iter]
  if (name == "gene_name" | name == "description") {
    corrected_names[iter]=name
  }
  else {
    corrected_names[iter]=desired[name]
  }
}
names(new_counts)=corrected_names

#Save file
write.table(new_counts,file=output_fname,row.names = F,sep="\t", quote=F)

#Summarize name revisions
num_changed = sum(corrected_names!=current_names)

#Revise renaming
df = data.frame("original"=current_names,"new"=colnames(new_counts))
to_remove = which(df[,1]==df[,2] & df[,1]%in%c("gene_name","description"))
df=df[-to_remove,]
same=(df$original==df$new)
notsame=!same
df$changed=notsame
row.names(df)=1:nrow(df)
fname = "results/sample_renaming_summary.txt"
write.table(df,file=fname,row.names = F,sep="\t", quote=F)
