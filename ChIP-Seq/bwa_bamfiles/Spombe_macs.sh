#!/bin/bash
################################################################################
# 
# 
# Created by David de la Cerda
# script name: Spombe_macs.sh
# 
# 
# 
# input: bam files
# output: xls files from macs2 with significant peaks and bed files with the closest genes 
# required software: R, bedtools, macs2
################################################################################
export PATH="$PATH:/home/ddelacer/.local/bin/"
source ~/.bashrc
#load R to get the pdf figure output
module load rhel7/R/4.0.2
export PATH=$PATH:/deac/peaseGrp/ddelacer/bedtools/bin/
source ~/.bashrc 
#min is 0.126e+8
#median is 0.134e+8
#max is 0.141e+8   
for sample in *bam
do
  echo $sample 
  describer=$(echo ${sample} |sed 's/.sorted.bam//')
  echo $describer
  macs2 callpeak -t $sample \
	-c /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq_11.16.2020/bwa_bamfiles/C1_wtbackground.sorted.bam \
 	-f BAMPE -g .141e+8 \
	-n ${describer} \
	--outdir /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq_11.16.2020/bwa_bamfiles/macs2_output 2> /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq_11.16.2020/bwa_bamfiles/macs2_output/${describer}--macs2.log
done 
cd /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq_11.16.2020/bwa_bamfiles/macs2_output
for sample in *_peaks.narrowPeak
do
    echo $sample
    describer=$(echo ${sample} |sed 's/_peaks.narrowPeak//')
    echo $describer
    closestBed -a $sample -b /deac/peaseGrp/ddelacer/pombe/data_collection/Schizosaccharomyces_pombe.ASM294v2.43.gtf > ${describer}_annotated.bed
done 
#tweaking file to get genes associated with peak of interest
for sample in *_annotated.bed
do
    echo $sample
    describer=$(echo ${sample} |sed 's/_annotated.bed//')
    #after merging it to the gtf file, sort by geene
    #keep columns 1 thru 19
    #parse out so that the final column has the gene name attached for the associated peak
    awk '$13 == "gene"' $sample | awk -v b=1 -v e=19 'BEGIN{FS=OFS="\t"} {for (i=b;i<=e;i++) printf "%s%s", $i, (i<e ? OFS : ORS)}' | awk -F "; " '$4 ~ "protein_coding" {print $1}' | sed 's/"//g' > ${describer}_peakgenes.bed
done
