#!/bin/bash
################################################################################
# 
# 
# Created by David de la Cerda
# script name: trim_adaptors.sh
# 
# 
# 
# input: untrimmed fq files and TruSeq2 adaptor fasta file
# output: trimmed fastq files used for read mapping
# required software: Trimmomatic-0.39
################################################################################
for file in *; do
    if [ -d "$file" ]; then
        # Will not run if no directories are available
        echo "$file"
        cd "$file"
        for sample in *_1.fq.gz
        do
            #echo $sample
            describer1=$(echo ${sample} |sed 's/.fq.gz//')
            echo $describer1
            describer2=$(echo ${sample} |sed 's/.fq.gz//' |sed 's/_1/_2/g')
            echo $describer2
            java -jar /deac/peaseGrp/ddelacer/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 ${describer1}.fq.gz ${describer2}.fq.gz ${describer1}.trim.fq.gz ${describer1}.untrim.fq.gz ${describer2}.trim.fq.gz ${describer2}.untrim.fq.gz SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:/deac/peaseGrp/ddelacer/pombe/pombe_ChIP/TruSeq2-PE.fa:2:40:15
            mv ${describer1}.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq
            mv ${describer2}.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq
            done
        cd ..
    fi
done