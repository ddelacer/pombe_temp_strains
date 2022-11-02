#!/bin/bash
################################################################################
# 
# 
# Created by David de la Cerda
# script name: fastQC.sh
# 
# 
# 
# input: fq files
# output: zip and html files
# required software: fastQC
################################################################################
for file in *; do
    if [ -d "$file" ]; then
        # Will not run if no directories are available
        echo "$file"
        cd "$file"
        for sample in *gz
        do
            echo $sample
            #describer=$(echo ${sample} |sed 's/.fastq.gz//')
            #echo $describer
            /deac/peaseGrp/ddelacer/FastQC/fastqc $sample -q -t 10 -d /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/fastqc_reports_11.16.2020
            mv *html /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/fastqc_reports_11.16.2020
            mv *zip /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/fastqc_reports_11.16.2020
            done
        cd ..
    fi
done