#!/bin/bash
################################################################################
# Data preparation S. pombe mutant samples published at: [DOI HERE]
# 
# Created by David de la Cerda
# script name: fastQC.sh
# 
# 
# 
# input:  trimmed fq.gz files
# output: FastQC reports
# required software: FastQC_v0.11.9
################################################################################
cd /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq
for sample in *gz
do
    echo $sample
    /deac/peaseGrp/ddelacer/FastQC/fastqc $sample -q -t 10 -d /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/report_fastqc_temp
    mv *html /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/report_fastqc_temp
    mv *zip /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/report_fastqc_temp
done