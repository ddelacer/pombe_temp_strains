#!/bin/bash
################################################################################
# Data preparation S. pombe mutant samples published at: [DOI HERE]
# 
# Created by David de la Cerda
# script name: bam_coverage_cpm.sh
# 
# 
# 
# input: bam files
# output: bigwig for visualization in genome viewer
# required software: python_3.7 deepTools_3.2.1
################################################################################
#CPM normalization 
for sample in *bam
do
    echo $sample
    describer=$(echo ${sample} |sed 's/Aligned.sortedByCoord.out.bam//')
    echo $describer
    bamCoverage -b $sample -bs 1 -p 10 --normalizeUsing CPM -o ${describer}_cpm.bw
    bamCoverage -b $sample -bs 1 -p 10 --normalizeUsing CPM --filterRNAstrand forward -o ${describer}_cpm_fwd.bw
    bamCoverage -b $sample -bs 1 -p 10 --normalizeUsing CPM --filterRNAstrand reverse -o ${describer}_cpm_rev.bw
done
mv *_cpm_fwd.bw forward/
mv *_cpm_rev.bw reverse/