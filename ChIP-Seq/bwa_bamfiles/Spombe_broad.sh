#!/bin/bash
################################################################################
# 
# 
# Created by David de la Cerda
# script name: Spombe_broad.sh
# 
# 
# 
# input: bam files
# output: xls files from macs2 with significant broad peaks and bed files with the closest genes 
# required software: bedtools, macs2
################################################################################
export PATH="$PATH:/home/ddelacer/.local/bin/"
source ~/.bashrc
export PATH=$PATH:/deac/peaseGrp/ddelacer/bedtools/bin/
source ~/.bashrc
for sample in *bam
do
  echo $sample 
  describer=$(echo ${sample} |sed 's/.sorted.bam//')
  echo $describer
  macs2 callpeak -t $sample \
	-c /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq_11.16.2020/bwa_bamfiles/C1_wtbackground.sorted.bam \
 	-f BAMPE -g .141e+8 \
	-n ${describer} \
    --broad \
	--outdir /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq_11.16.2020/bwa_bamfiles/macs2_broad_output 2> /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq_11.16.2020/bwa_bamfiles/macs2_broad_output/${describer}--macs2.log
done 
cd /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq_11.16.2020/bwa_bamfiles/macs2_broad_output
#annotate the broad peaks file
for sample in *_peaks.broadPeak
do
    echo $sample
    describer=$(echo ${sample} |sed 's/_peaks.broadPeak//')
    echo $describer
    closestBed -a $sample -b /deac/peaseGrp/ddelacer/pombe/data_collection/Schizosaccharomyces_pombe.ASM294v2.43.gtf > ${describer}_broad_annotated.bed
done 
for sample in *_broad_annotated.bed
do
    echo $sample
    describer=$(echo ${sample} |sed 's/_broad_annotated.bed//')
#    #after merging it to the gtf file, sort by geene
#    #keep columns 1 thru 19
#    #parse out so that the final column has the gene name attached for the associated peak
    awk '$13 == "gene"' $sample | awk -v b=1 -v e=19 'BEGIN{FS=OFS="\t"} {for (i=b;i<=e;i++) printf "%s%s", $i, (i<e ? OFS : ORS)}' | awk -F "; " '$4 ~ "protein_coding" {print $1}' | sed 's/"//g' > ${describer}_broadpeakgenes.bed
done
