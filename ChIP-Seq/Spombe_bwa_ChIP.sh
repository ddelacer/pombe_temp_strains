#!/bin/bash
################################################################################
# 
# 
# Created by David de la Cerda
# script name: Spombe_bwa_ChIP.sh
# 
# 
# 
# input: fastq files
# output: bam files for analysis, fasta to map to is in ChIP-Seq directory
# required software: bwa_0.7.17 samtools_1.10
################################################################################
#module load rhel7/bwa/0.7.17
#module load rhel7/gcc/10.1.0 rhel7/samtools/1.10
A1_lsd2ftp="/deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/A1_1.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/A1_2.trim.fq.gz"
A2_lsd2deltacftp="/deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/A2_1.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/A2_2.trim.fq.gz"
A3_lsd2ftpclr4del="/deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/A3_1.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/A3_2.trim.fq.gz"
A4_lsd2ftpset1del="/deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/A4_1.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/A4_2.trim.fq.gz"
A5_lsd2ftpclr4set1del="/deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/A5_1.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/A5_2.trim.fq.gz"
C1_wtbackground="/deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/C1_1.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/C1_2.trim.fq.gz"
C2_lsd1ftp="/deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/C2_1.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/C2_2.trim.fq.gz"
C3_lsd1deltahmgftp="/deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/C3_1.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/C3_2.trim.fq.gz"
C4_lsd1ftpclr4del="/deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/C4_1.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/C4_2.trim.fq.gz"
C5_lsd1ftpset1del="/deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/C5_1.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/C5_2.trim.fq.gz"
C6_lsd1ftpclr4set1del="/deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/C6_1.trim.fq.gz /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq/C6_2.trim.fq.gz"
bwa index /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes.fa
bwa mem  -t 1 /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes.fa $A1_lsd2ftp > A1_lsd2ftp.sam 
bwa mem  -t 1 /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes.fa $A2_lsd2deltacftp > A2_lsd2deltacftp.sam
bwa mem  -t 1 /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes.fa $A3_lsd2ftpclr4del > A3_lsd2ftpclr4del.sam
bwa mem  -t 1 /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes.fa $A4_lsd2ftpset1del > A4_lsd2ftpset1del.sam
bwa mem  -t 1 /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes.fa $A5_lsd2ftpclr4set1del > A5_lsd2ftpclr4set1del.sam
bwa mem  -t 1 /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes.fa $C1_wtbackground > C1_wtbackground.sam
bwa mem  -t 1 /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes.fa $C2_lsd1ftp > C2_lsd1ftp.sam
bwa mem  -t 1 /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes.fa $C3_lsd1deltahmgftp > C3_lsd1deltahmgftp.sam
bwa mem  -t 1 /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes.fa $C4_lsd1ftpclr4del > C4_lsd1ftpclr4del.sam
bwa mem  -t 1 /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes.fa $C5_lsd1ftpset1del > C5_lsd1ftpset1del.sam
bwa mem  -t 1 /deac/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes.fa $C6_lsd1ftpclr4set1del > C6_lsd1ftpclr4set1del.sam
for sample in *.sam
do
    echo $sample
    describer=$(echo ${sample} |sed 's/.sam//')
    samtools sort $sample -o ${describer}.sorted.bam 
done
#index new sorted bam files
for bamfile in *.sorted.bam
do
    echo $bamfile
    samtools index $bamfile
done
mv *bam bwa_bam_files/
cd /bwa_bam_files/
for sample in *.sorted.bam
do
    echo $sample
    describer=$(echo ${sample} |sed 's/.sorted.bam//')
    echo $describer
    samtools stats $sample > ${describer}.stats.out
done
#remove all the unnecessary and large sam files
rm *.sam
for sample in *.stats.out
do 
    echo $sample
    describer=$(echo ${sample} |sed 's/.stats.out//')
    echo $describer
    grep ^SN $sample | cut -f 2- > ${describer}.summary.out
done 
#combine all the summary files for comparison
paste *.summary.out > all.summary.txt
mv all.summary.txt ..