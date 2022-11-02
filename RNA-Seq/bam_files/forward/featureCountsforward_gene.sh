#!/bin/bash
################################################################################
# Data preparation S. pombe mutant samples published at: [DOI HERE]
# 
# Created by David de la Cerda
# script name: featureCounts_gene.sh
# 
# 
# 
# input:  bam files
# output: gene counts file for statistical analysis in directory called gene_counts
# required software: featureCounts_1.6.3
################################################################################
#making count tables for forward strand from RNA-Seq
for sample in *bam
do
    echo $sample
    describer=$(echo ${sample} |sed 's/_forward.bam//')
    echo $describer
    #gene
    /deac/bio/peaseGrp/ddelacer/subread-1.6.3-source/bin/featureCounts -t gene -g gene_id -a /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf -s 2 -G /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes_.fa -O -M --fraction -T 8 -p -o ${describer}_gene.txt $sample
    cut -f1,7-8 ${describer}_gene.txt | sed -e '1d' - > ${describer}_gene_final.txt
done
mv *gene* gene_counts/