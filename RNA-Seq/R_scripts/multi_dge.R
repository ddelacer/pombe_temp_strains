################################################################################
# RNA-Seq analysis for various S. pombe samples published at: [DOI]
# Compraisons consist of:
# temp mutants 37C:  lsd1/2 set1
# nontemp mutants: lsd1/2 clr4 clr4/lsd1 set1 set1/lsd2 clr6 clr6/lsd1 clr6/lsd2 sir2 sir2/lsd1 sir2/lsd2 
# wt temp to nontemp
# Created by David de la Cerda
# script name: multi_dge.R
# 
# 
# 
# input: gene quantification file generated from  featureCounts_gene.sh, featureCountsforward_gene.sh, featureCountsreverse_gene.sh
# output:  DGE output files, Top gene output files for GO analysis
# required software: R version 4.0.4
################################################################################

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
gc()

library(edgeR)
library(xlsx)

###Functions
#read in all genes featureCounts table
readit = function(filename,dfout){
  dfout = read.table(file=filename,header=T,row.names=1,sep="\t")
  names(dfout) = sub("Aligned.sortedByCoord.out.bam","",names(dfout))
  dfout = dfout[,!grepl("^Geneid",names(dfout))]
  dfout = floor(dfout)
  return(dfout)
}

#read in forward featureCounts table
readitfwd = function(filename,dfout){
  dfout = read.table(file=filename,header=T,row.names=1,sep="\t")
  names(dfout) = sub("_forward.bam","",names(dfout))
  dfout = dfout[,!grepl("^Geneid",names(dfout))]
  dfout = floor(dfout)
  return(dfout)
}

#read in reverse featureCounts table
readitrev = function(filename,dfout){
  dfout = read.table(file=filename,header=T,row.names=1,sep="\t")
  names(dfout) = sub("_reverse.bam","",names(dfout))
  dfout = dfout[,!grepl("^Geneid",names(dfout))]
  dfout = floor(dfout)
  return(dfout)
}

#DGE analysis using edgeR
dgplus = function(df,group1,group2,grouping){
  grp1 = group1
  grp2 = group2
  #subset original dataframe by merging string groups
  temp1 = subset(df,select=grp1)
  temp2 = subset(df,select=grp2)
  #new subset dataframe below
  df = transform(merge(temp1,temp2,by="row.names"),row.names=Row.names,Row.names=NULL)
  
  y <- df
  y_full <- DGEList(counts=y,group = grouping)
  keep <- rowSums(cpm(y_full) > 0) >= 2
  lm1.y <- y_full[keep, keep.lib.sizes =FALSE]
  lm1.y <- calcNormFactors(lm1.y, method="TMM")
  lm1.tmm = as.data.frame(lm1.y$counts)
  lm1.design <- model.matrix( ~ as.numeric(group==1), data = lm1.y$samples)
  lm1.v <- voom(lm1.y, lm1.design, plot=FALSE)
  lm1.fit <- lmFit(lm1.v, lm1.design)
  lm1.fit.c <- eBayes(lm1.fit)
  lm1.top <- topTable(lm1.fit.c, sort="P",number='all')
  lm1.top.merge = transform(merge(lm1.top,df,by="row.names"),row.names=Row.names,Row.names=NULL)
  lm1.top.merge = lm1.top.merge[order(lm1.top.merge$adj.P.Val),]
  return(lm1.top.merge)
}

#This function outputs excel file from generated list
list2excelout = function(datalist,sheetnamelist,filename){
  wb <- createWorkbook()
  datas <- datalist
  sheetnames = names(sheetnamelist)
  sheets <- lapply(sheetnames, createSheet, wb = wb)
  void <- Map(addDataFrame, datas, sheets)
  saveWorkbook(wb, file = filename)
}

#subsetting significant genes into gene lists for GO
#these files can be input into Ontologizer (Java)
#Right now significance is hard-coded to specific P-Value threshold
list2subsettable = function(datalist,sheetnamelist,fileextension){
  #first element of datalist and sheetname list are empty, that's why -1
  sublist = lapply(datalist[-1], function(x) {subset(x,adj.P.Val < .0001)} )
  for (strain in seq_along(sublist)){
    filename = paste0(names(sheetnamelist[-1])[strain], fileextension)
    write.table( rownames(sublist[[strain]]) ,filename,quote=F,sep="\t",row.names=F,col.names=F)
  }
  
}

list2subsettable_custom = function(datalist,sheetnamelist,fileextension){
  #first element of datalist and sheetname list are empty, that's why -1
  sublist = lapply(datalist, function(x) {subset(x,adj.P.Val < .0001)} )
  for (strain in seq_along(sublist)){
    filename = paste0(names(sheetnamelist)[strain], fileextension)
    write.table( rownames(sublist[[strain]]) ,filename,quote=F,sep="\t",row.names=F,col.names=F)
  }
  
}

#Output the population list of genes for each contrast tested
list2poptable = function(datalist,sheetnamelist,fileextension){
  #first element of datalist and sheetname list are empty, that's why -1
  #sublist = lapply(datalist[-1], function(x) {subset(x,adj.P.Val < .0001)} )
  sublist=datalist[-1]
  for (strain in seq_along(sublist)){
    filename = paste0(names(sheetnamelist[-1])[strain], fileextension)
    write.table( rownames(sublist[[strain]]) ,filename,quote=F,sep="\t",row.names=F,col.names=F)
  }
  
}

list2poptable_custom = function(datalist,sheetnamelist,fileextension){
  #first element of datalist and sheetname list are empty, that's why -1
  #sublist = lapply(datalist[-1], function(x) {subset(x,adj.P.Val < .0001)} )
  sublist=datalist
  for (strain in seq_along(sublist)){
    filename = paste0(names(sheetnamelist)[strain], fileextension)
    write.table( rownames(sublist[[strain]]) ,filename,quote=F,sep="\t",row.names=F,col.names=F)
  }
  
}


#Gene quantification tables from featureCounts
pombe_geneOfrac_all = readit("Spombe_gene_withO_fraction_all.txt")
pombe_geneOfrac_fwd = readitfwd("Spombe_gene_withO_fraction_fwd.txt")
pombe_geneOfrac_rev = readitrev("Spombe_gene_withO_fraction_rev.txt")

nontemp_lst = list(
  "wt" = "A1_WT_1|B1_WT_2",
  "lsd1" = "A3|B3",
  "lsd2" = "A5|B5",
  "clr4" = "A7|B7",
  "clr4_lsd1" = "A8|B8",
  "set1" = "A9|B9",
  "set1_lsd2" = "A11|B11",
  "clr6" = "A12|B12",
  "clr6_lsd1" = "A13|B13",
  "clr6_lsd2" = "A14|B14",
  "sir2" = "A15|B15",
  "sir2_lsd1" = "A16|B16",
  "sir2_lsd2" = "A17|B17"
  
)

temp_lst = list(
  "wt_37" = "A2|B2",
  "lsd1_37" = "A4|B4",
  "lsd2_37" = "A6|B6",
  "set1_37" = "A10|B10"
)

custom_lst = list(
  "wt" = "A1_WT_1|B1_WT_2", #30 degrees
  "wt_37" = "A2|B2",
  "lsd1" = "A3|B3",
  "lsd1_37" = "A4|B4",
  "lsd2" = "A5|B5",
  "lsd2_37" = "A6|B6",
  "set1" = "A9|B9",
  "set1_37" = "A10|B10"
)

#custom output
cusdgeoutput = vector("list",(length(custom_lst)/2) )
#even numbers refer to mutant, odd to control using this scheme
for (strain in (seq(2,(length(custom_lst)),by=2 ) ) ){
  controls = grep(custom_lst[[strain-1]],names(pombe_geneOfrac_all),value=T)
  mutants = grep(custom_lst[[strain]],names(pombe_geneOfrac_all),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  cusdgeoutput[[strain]]=dgplus(pombe_geneOfrac_all,controls,mutants,compare_grp)
  
}

customnames = list(
  "wt_30_wt_37" = "wt_30/wt_37" ,
  "lsd1_30_lsd1_37" = "lsd1_30/lsd1_37",
  "lsd2_30_lsd2_37" ="lsd2_30/lsd2_37",
  "set1_30_set1_37" = "set1_30/set1_37"
)
#export the output
list2excelout(cusdgeoutput,custom_lst,"customtemp_rnaseq.xlsx") 
#remove null values generated accidentally
cusdgeoutput[sapply(cusdgeoutput, is.null)] <- NULL
#significant gene list
list2subsettable_custom(cusdgeoutput,customnames,"_customtemp_all_list.txt")
#pop gene list for GO
list2poptable_custom(cusdgeoutput,customnames,"_customtemp_all_pop.txt")


#all genes
nontempdgeout = vector("list",(length(nontemp_lst)-1) )
for (strain in (seq(2:(length(nontemp_lst)) )+1) ){
  controls = grep(nontemp_lst[[1]],names(pombe_geneOfrac_all),value=T)
  mutants = grep(nontemp_lst[[strain]],names(pombe_geneOfrac_all),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  nontempdgeout[[strain]]=dgplus(pombe_geneOfrac_all,controls,mutants,compare_grp)
  
}

#forward strand
nontempdgeout_fwd = vector("list",(length(nontemp_lst)-1) )
for (strain in (seq(2:(length(nontemp_lst)) )+1) ){
  controls = grep(nontemp_lst[[1]],names(pombe_geneOfrac_fwd),value=T)
  mutants = grep(nontemp_lst[[strain]],names(pombe_geneOfrac_fwd),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  nontempdgeout_fwd[[strain]]=dgplus(pombe_geneOfrac_fwd,controls,mutants,compare_grp)
  
}

#reverse strand
nontempdgeout_rev = vector("list",(length(nontemp_lst)-1) )
for (strain in (seq(2:(length(nontemp_lst)) )+1) ){
  controls = grep(nontemp_lst[[1]],names(pombe_geneOfrac_rev),value=T)
  mutants = grep(nontemp_lst[[strain]],names(pombe_geneOfrac_rev),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  nontempdgeout_rev[[strain]]=dgplus(pombe_geneOfrac_rev,controls,mutants,compare_grp)
  
}

#Export into Excel files
wb <- createWorkbook()
datas <- nontempdgeout
sheetnames = names(nontemp_lst)
sheets <- lapply(sheetnames, createSheet, wb = wb)
void <- Map(addDataFrame, datas, sheets)
saveWorkbook(wb, file = "nontempall_rnaseq.xlsx")

#list2excelout(nontempdgeout_fwd,nontemp_lst,"nontempfwd_rnaseq.xlsx")  
#list2excelout(nontempdgeout_rev,nontemp_lst,"nontemprev_rnaseq.xlsx")  

list2subsettable(nontempdgeout,nontemp_lst,"_nontemp_all_list.txt")
list2subsettable(nontempdgeout_fwd,nontemp_lst,"_nontemp_fwd_list.txt")
list2subsettable(nontempdgeout_rev,nontemp_lst,"_nontemp_rev_list.txt")

#temperature mutants
tempdgeout = vector("list",(length(temp_lst)-1) )
for (strain in (seq(2:(length(temp_lst)) )+1) ){
  controls = grep(temp_lst[[1]],names(pombe_geneOfrac_all),value=T)
  mutants = grep(temp_lst[[strain]],names(pombe_geneOfrac_all),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  tempdgeout[[strain]]=dgplus(pombe_geneOfrac_all,controls,mutants,compare_grp)
  
}

#forward strand
tempdgeout_fwd = vector("list",(length(temp_lst)-1) )
for (strain in (seq(2:(length(temp_lst)) )+1) ){
  controls = grep(temp_lst[[1]],names(pombe_geneOfrac_fwd),value=T)
  mutants = grep(temp_lst[[strain]],names(pombe_geneOfrac_fwd),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  tempdgeout_fwd[[strain]]=dgplus(pombe_geneOfrac_fwd,controls,mutants,compare_grp)
  
}

#reverse strand
tempdgeout_rev = vector("list",(length(temp_lst)-1) )
for (strain in (seq(2:(length(temp_lst)) )+1) ){
  controls = grep(temp_lst[[1]],names(pombe_geneOfrac_rev),value=T)
  mutants = grep(temp_lst[[strain]],names(pombe_geneOfrac_rev),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  tempdgeout_rev[[strain]]=dgplus(pombe_geneOfrac_rev,controls,mutants,compare_grp)
  
}

#output here
wb <- createWorkbook()
datas <- tempdgeout
sheetnames = names(temp_lst)
sheets <- lapply(sheetnames, createSheet, wb = wb)
void <- Map(addDataFrame, datas, sheets)
saveWorkbook(wb, file = "temp_rnaseq.xlsx")

#list2excelout(tempdgeout_fwd,temp_lst,"tempfwd_rnaseq.xlsx")
#list2excelout(tempdgeout_rev,temp_lst,"temprev_rnaseq.xlsx")

list2subsettable(tempdgeout,temp_lst,"_temp_all_list.txt")
list2subsettable(tempdgeout_fwd,temp_lst,"_temp_fwd_list.txt")
list2subsettable(tempdgeout_rev,temp_lst,"_temp_rev_list.txt")

#create population files for GO
#each tested set has a different number of testable genes, might affect GO output
list2poptable(nontempdgeout,nontemp_lst,"_nontemp_all_pop.txt")
list2poptable(nontempdgeout_fwd,nontemp_lst,"_nontemp_fwd_pop.txt")
list2poptable(nontempdgeout_rev,nontemp_lst,"_nontemp_rev_pop.txt")
list2poptable(tempdgeout,temp_lst,"_temp_all_pop.txt")
list2poptable(tempdgeout_fwd,temp_lst,"_temp_fwd_pop.txt")
list2poptable(tempdgeout_rev,temp_lst,"_temp_rev_pop.txt")

#compare wt nontemp vs temp controls
wtdgeout = vector("list",(length(1)) )
for (strain in (seq(1:1)) ){
  mutants = grep(temp_lst[[1]],names(pombe_geneOfrac_all),value=T)
  controls = grep(nontemp_lst[[1]],names(pombe_geneOfrac_all),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  wtdgeout[[strain]]=dgplus(pombe_geneOfrac_all,controls,mutants,compare_grp)
  
}

#forward strand
wtdgeout_fwd = vector("list",(length(1)) )
for (strain in (seq(1:1)) ){
  mutants = grep(temp_lst[[1]],names(pombe_geneOfrac_fwd),value=T)
  controls = grep(nontemp_lst[[1]],names(pombe_geneOfrac_fwd),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  wtdgeout_fwd[[strain]]=dgplus(pombe_geneOfrac_fwd,controls,mutants,compare_grp)
  
}

#reverse strand
wtdgeout_rev = vector("list",(length(1)) )
for (strain in (seq(1:1)) ){
  mutants = grep(temp_lst[[1]],names(pombe_geneOfrac_rev),value=T)
  controls = grep(nontemp_lst[[1]],names(pombe_geneOfrac_rev),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  wtdgeout_rev[[strain]]=dgplus(pombe_geneOfrac_rev,controls,mutants,compare_grp)
  
}

#can combine wt output into one excel file for simplicity
wtdgeout_combined = c(wtdgeout,wtdgeout_fwd,wtdgeout_rev)

#output here
wb <- createWorkbook()
datas <- wtdgeout_combined
sheetnames = c("wt_all","wt_forward","wt_reverse")
sheets <- lapply(sheetnames, createSheet, wb = wb)
void <- Map(addDataFrame, datas, sheets)
saveWorkbook(wb, file = "wt_dge_rnaseq.xlsx")
