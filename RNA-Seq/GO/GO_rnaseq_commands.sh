#!/bin/bash
################################################################################
# GO analysis of S. pombe mutant samples published at: [DOI HERE]
# 
# Created by David de la Cerda
# script name: GO_rnaseq_commands.sh
# 
# 
# 
# input: Output DGE lists and pop files from multi_dge.R
# output: GO tables and annotation files
# required software: Ontologizer_2.1
#################################################################################
These commands should perform GO on the given list of DEGs per comparison
for sample in *_list.txt
#for sample in *_customtemp_all_list.txt
do
    echo $sample
    describer=$(echo ${sample} |sed 's/_list.txt//')
    echo $describer
    java -jar /deac/bio/peaseGrp/ddelacer/Ontologizer.jar -a pombase.gaf -g go.obo -s $sample -p ${describer}_pop.txt -c Parent-Child-Union -n -m Benjamini-Hochberg
done