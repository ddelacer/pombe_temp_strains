#!/bin/bash
################################################################################
# Data preparation S. pombe mutant samples published at: [DOI HERE]
# 
# Created by David de la Cerda
# script name: bam_index.sh
# 
# 
# 
# input:  bam files
# output: bam index files
# required software: samtools_1.10
################################################################################
module load rhel7/samtools/1.10
for sample in *bam
do
    echo $sample
    samtools index $sample
done