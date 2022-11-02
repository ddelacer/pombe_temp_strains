#!/bin/bash
################################################################################
# Data preparation S. pombe mutant samples published at: [DOI HERE]
# 
# Created by David de la Cerda
# script name: Spombe_strand_neg.sh
# 
# 
# 
# input:  bam files
# output: forward strand bam files in directory called reverse
# required software: python_3.7 samtools_1.10
################################################################################
module load rhel7/compilers/intel-2018-lp64
module load rhel7/samtools/1.10
module load rhel7/python/3.7.0
export PATH="$PATH:/home/ddelacer/.local/bin/"
source ~/.bashrc
#http://broadinstitute.github.io/picard/explain-flags.html
#using script adapted from:
#https://www.biostars.org/p/92935/
for sample in *bam
#negative strand
do
    describer=$(echo ${sample} |sed 's/Aligned.sortedByCoord.out.bam//')
    echo $describer
    samtools view -b -f 144 -o ${describer}_rvs1.bam $sample
    samtools view -b -f 64 -F 16 -o ${describer}_rvs2.bam $sample
    
done
mv *_rvs?.bam reverse/
cd reverse/
for sample in *bam
do
    samtools index $sample
done
#merge the reverse options together
for sample in *_rvs1.bam
do
    describer=$(echo ${sample} |sed 's/_rvs1.bam//')
    echo $describer
    samtools merge -f ${describer}_reverse.bam ${describer}_rvs1.bam ${describer}_rvs2.bam 
done
for sample in *reverse.bam
do
    samtools index $sample
done
rm *rvs*