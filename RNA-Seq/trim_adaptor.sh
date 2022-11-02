#!/bin/bash
################################################################################
# Data preparation S. pombe mutant samples published at: [DOI HERE]
# 
# Created by David de la Cerda
# script name: trim_adaptor.sh
# 
# 
# 
# input: fq.gz files uploaded to GEO [GSE HERE]
# output: trimmed fq.gz files
# required software: trimmomatic_0.39
################################################################################
module load rhel7/compilers/intel-2018-lp64
export PATH="$PATH:/home/ddelacer/.local/bin/"
source ~/.bashrc
#finalized version of trimming script
#files were stored in individual directories originally
#can adapt script to work if all files are within one directory if 
#downloaded from GEO 
for file in *; do
    if [ -d "$file" ]; then
        # Will not run if no directories are available
        echo "$file"
        cd "$file"
        for sample in *_1.fq.gz
        do
            describer1=$(echo ${sample} |sed 's/.fq.gz//')
            echo $describer1
            describer2=$(echo ${sample} |sed 's/.fq.gz//' |sed 's/_1/_2/g')
            echo $describer2
            java -jar /deac/peaseGrp/ddelacer/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 ${describer1}.fq.gz ${describer2}.fq.gz ${describer1}.trim.fq.gz ${describer1}.untrim.fq.gz ${describer2}.trim.fq.gz ${describer2}.untrim.fq.gz SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/raw_data/TruSeq2-PE.fa:2:40:15
            mv ${describer1}.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq
            mv ${describer2}.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq
            done
        cd ..
    fi
done