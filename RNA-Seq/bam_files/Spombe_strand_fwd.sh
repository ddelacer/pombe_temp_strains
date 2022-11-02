#!/bin/bash
################################################################################
# Data preparation S. pombe mutant samples published at: [DOI HERE]
# 
# Created by David de la Cerda
# script name: Spombe_strand_fwd.sh
# 
# 
# 
# input:  bam files
# output: forward strand bam files in directory called forward
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
#forward strand
do
    describer=$(echo ${sample} |sed 's/Aligned.sortedByCoord.out.bam//')
    echo $describer
    #samtools view -f 17 $sample -o ${describer}_forward.bam
    samtools view -b -f 128 -F 16 -o ${describer}_fwd1.bam $sample
    samtools view -b -f 80 -o ${describer}_fwd2.bam $sample
    
done
mv *_fwd?.bam forward/
cd forward/
for sample in *bam
do
    samtools index $sample
done
#merge the forward options together
for sample in *_fwd1.bam
do
    describer=$(echo ${sample} |sed 's/_fwd1.bam//')
    echo $describer
    samtools merge -f ${describer}_forward.bam ${describer}_fwd1.bam ${describer}_fwd2.bam 
done
for sample in *forward.bam
do
    samtools index $sample
done
rm *fwd*
