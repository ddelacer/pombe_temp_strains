#!/bin/bash
################################################################################
# Data preparation S. pombe mutant samples published at: [DOI HERE]
# 
# Created by David de la Cerda
# script name: Spombe_map_temp.sh
# 
# 
# 
# input: trimmed fq.gz files
# output: mapped bam files for analysis
# required software: star_2.7.2d
################################################################################
module load rhel7/compilers/intel-2018-lp64
cd/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/bam_files
#Give all the file names variable names and then do mapping, each variable has file for each paired-end read
A1_WT_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A1_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A1_2.trim.fq.gz"
A2_WT_37C_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A2_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A2_2.trim.fq.gz"
A3_lsd1_dHMG_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A3_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A3_2.trim.fq.gz"
A4_lsd1_dHMG_37C_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A4_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A4_2.trim.fq.gz"
A5_lsd2_dC_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A5_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A5_2.trim.fq.gz"
A6_lsd2_dC_37C_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A6_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A6_2.trim.fq.gz"
A7_clr4d_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A7_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A7_2.trim.fq.gz"
A8_clr4d_lsd1_dHMG_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A8_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A8_2.trim.fq.gz"
A9_set1d_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A9_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A9_2.trim.fq.gz"
A10_set1d_37C_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A10_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A10_2.trim.fq.gz"
A11_set1d_lsd2_dC_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A11_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A11_2.trim.fq.gz"
A12_clr6_1_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A12_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A12_2.trim.fq.gz"
A13_clr6_1_lsd1_dHMG_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A13_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A13_2.trim.fq.gz"
A14_clr6_1_lsd2_dC_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A14_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A14_2.trim.fq.gz"
A15_sir2d_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A15_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A15_2.trim.fq.gz"
A16_sir2d_lsd1_dHMG_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A16_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A16_2.trim.fq.gz"
A17_sir2d_lsd2_dC_1="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A17_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/A17_2.trim.fq.gz"
B1_WT_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B1_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B1_2.trim.fq.gz"
B2_WT_37C_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B2_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B2_2.trim.fq.gz"
B3_lsd1_dHMG_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B3_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B3_2.trim.fq.gz"
B4_lsd1_dHMG_37C_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B4_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B4_2.trim.fq.gz"
B5_lsd2_dC_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B5_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B5_2.trim.fq.gz"
B6_lsd2_dC_37C_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B6_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B6_2.trim.fq.gz"
B7_clr4d_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B7_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B7_2.trim.fq.gz"
B8_clr4d_lsd1_dHMG_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B8_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B8_2.trim.fq.gz"
B9_set1d_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B9_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B9_2.trim.fq.gz"
B10_set1d_37C_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B10_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B10_2.trim.fq.gz"
B11_set1d_lsd2_dC_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B11_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B11_2.trim.fq.gz"
B12_clr6_1_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B12_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B12_2.trim.fq.gz"
B13_clr6_1_lsd1_dHMG_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B13_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B13_2.trim.fq.gz"
B14_clr6_1_lsd2_dC_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B14_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B14_2.trim.fq.gz"
B15_sir2d_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B15_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B15_2.trim.fq.gz"
B16_sir2d_lsd1_dHMG_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B16_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B16_2.trim.fq.gz"
B17_sir2d_lsd2_dC_2="/deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B17_1.trim.fq.gz /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/trimmed_fq/B17_2.trim.fq.gz"
#Genome Generate step
module load rhel7/star/2.7.2d
STAR --runMode genomeGenerate --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --genomeFastaFiles /deac/bio/peaseGrp/ddelacer/pombe_genome/pombe-genome-starIndex/Schizosaccharomyces_pombe_all_chromosomes_.fa --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --sjdbOverhang 149 --genomeSAindexNbases 10
echo $(basename $A1_WT_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A1_WT_1 --outFileNamePrefix A1_WT_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A2_WT_37C_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A2_WT_37C_1 --outFileNamePrefix A2_WT_37C_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A3_lsd1_dHMG_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A3_lsd1_dHMG_1 --outFileNamePrefix A3_lsd1_dHMG_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A4_lsd1_dHMG_37C_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A4_lsd1_dHMG_37C_1 --outFileNamePrefix A4_lsd1_dHMG_37C_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A5_lsd2_dC_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A5_lsd2_dC_1 --outFileNamePrefix A5_lsd2_dC_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A6_lsd2_dC_37C_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A6_lsd2_dC_37C_1 --outFileNamePrefix A6_lsd2_dC_37C_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A7_clr4d_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A7_clr4d_1 --outFileNamePrefix A7_clr4d_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A8_clr4d_lsd1_dHMG_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A8_clr4d_lsd1_dHMG_1 --outFileNamePrefix A8_clr4d_lsd1_dHMG_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A9_set1d_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A9_set1d_1 --outFileNamePrefix A9_set1d_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A10_set1d_37C_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A10_set1d_37C_1 --outFileNamePrefix A10_set1d_37C_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A11_set1d_lsd2_dC_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A11_set1d_lsd2_dC_1 --outFileNamePrefix A11_set1d_lsd2_dC_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A12_clr6_1_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A12_clr6_1_1 --outFileNamePrefix A12_clr6_1_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A13_clr6_1_lsd1_dHMG_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A13_clr6_1_lsd1_dHMG_1 --outFileNamePrefix A13_clr6_1_lsd1_dHMG_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A14_clr6_1_lsd2_dC_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A14_clr6_1_lsd2_dC_1 --outFileNamePrefix A14_clr6_1_lsd2_dC_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A15_sir2d_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A15_sir2d_1 --outFileNamePrefix A15_sir2d_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A16_sir2d_lsd1_dHMG_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A16_sir2d_lsd1_dHMG_1 --outFileNamePrefix A16_sir2d_lsd1_dHMG_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $A17_sir2d_lsd2_dC_1)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $A17_sir2d_lsd2_dC_1 --outFileNamePrefix A17_sir2d_lsd2_dC_1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B1_WT_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B1_WT_2 --outFileNamePrefix B1_WT_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B2_WT_37C_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B2_WT_37C_2 --outFileNamePrefix B2_WT_37C_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B3_lsd1_dHMG_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B3_lsd1_dHMG_2 --outFileNamePrefix B3_lsd1_dHMG_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B4_lsd1_dHMG_37C_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B4_lsd1_dHMG_37C_2 --outFileNamePrefix B4_lsd1_dHMG_37C_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B5_lsd2_dC_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B5_lsd2_dC_2 --outFileNamePrefix B5_lsd2_dC_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B6_lsd2_dC_37C_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B6_lsd2_dC_37C_2 --outFileNamePrefix B6_lsd2_dC_37C_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B7_clr4d_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B7_clr4d_2 --outFileNamePrefix B7_clr4d_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B8_clr4d_lsd1_dHMG_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B8_clr4d_lsd1_dHMG_2 --outFileNamePrefix B8_clr4d_lsd1_dHMG_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B9_set1d_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B9_set1d_2 --outFileNamePrefix B9_set1d_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B10_set1d_37C_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B10_set1d_37C_2 --outFileNamePrefix B10_set1d_37C_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B11_set1d_lsd2_dC_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B11_set1d_lsd2_dC_2 --outFileNamePrefix B11_set1d_lsd2_dC_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B12_clr6_1_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B12_clr6_1_2 --outFileNamePrefix B12_clr6_1_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B13_clr6_1_lsd1_dHMG_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B13_clr6_1_lsd1_dHMG_2 --outFileNamePrefix B13_clr6_1_lsd1_dHMG_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B14_clr6_1_lsd2_dC_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B14_clr6_1_lsd2_dC_2 --outFileNamePrefix B14_clr6_1_lsd2_dC_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B15_sir2d_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B15_sir2d_2 --outFileNamePrefix B15_sir2d_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B16_sir2d_lsd1_dHMG_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B16_sir2d_lsd1_dHMG_2 --outFileNamePrefix B16_sir2d_lsd1_dHMG_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244
echo $(basename $B17_sir2d_lsd2_dC_2)
STAR --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com --readFilesCommand zcat --readFilesIn $B17_sir2d_lsd2_dC_2 --outFileNamePrefix B17_sir2d_lsd2_dC_2 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/Schizosaccharomyces_pombe.ASM294v2.43.gtf --runThreadN 4  --limitBAMsortRAM 7787218244