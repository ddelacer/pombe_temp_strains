#!/bin/bash
################################################################################
# Data preparation lsd1/lsd2 samples published at: https://doi.org/10.3390/cells9040955
# 
# Created by David de la Cerda
# script name: chip_newcompare_3.31.2022.sh
# 
# 
# 
# input: log2 bigwig files and narrow peak cluster bed files
# output: matrix text file and png heatmap image of peaks across samples 
# required software: python_3.7 deepTools_3.2.1
################################################################################
export PATH="$PATH:/home/ddelacer/.local/bin/"
source ~/.bashrc
#making bed file of genes from MACS2 output
cd /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq_11.16.2020/bwa_bamfiles/macs2_output/
for sample in *peakgenes.bed
do
    echo $sample
    describer=$(echo ${sample} |sed 's/_peakgenes.bed//')
    echo $describer
    #find the intersections between peak genes keeping original entries of both files
    intersectBed -a $sample -b /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/trimmed_fq_11.16.2020/bwa_bamfiles/macs2_output/${sample} > /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/combined_peaks/${describer}_narrowcombine.bed
    closestBed -a /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/combined_peaks/${describer}_narrowcombine.bed -b /deac/peaseGrp/ddelacer/pombe/data_collection/Schizosaccharomyces_pombe.ASM294v2.43.gtf > /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/combined_peaks/${describer}_annotated.bed 
done
mv *annotated_bed /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/combined_peaks/
cd /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/combined_peaks/
for sample in *_annotated.bed
do
    echo $sample
    describer=$(echo ${sample} |sed 's/_annotated.bed//')
    awk '$23 == "gene"' $sample | awk -v b=1 -v e=19 'BEGIN{FS=OFS="\t"} {for (i=b;i<=e;i++) printf "%s%s", $i, (i<e ? OFS : ORS)}' > ${describer}_finalcombedpeakgenes.bed
done
#need to make customized gtf files because deeptools has specific requirements
#make a list that is separated by semicolons
cd /deac/peaseGrp/ddelacer/pombe/pombe_ChIP/combined_peaks/
for sample in *_finalcombedpeakgenes.bed
do
    echo $sample
    describer=$(echo ${sample} |sed 's/_finalcombedpeakgenes.bed//')
    echo $describer
    awk '{print $(NF)}' $sample |sort | uniq > ${describer}_genelist.txt
    grep -f ${describer}_genelist.txt /deac/peaseGrp/ddelacer/pombe/data_collection/Schizosaccharomyces_pombe.ASM294v2.43.gtf > ${describer}.gtf
    #add extra line because deeptools always errors line one for some reason
    head -n 1 ${describer}.gtf | cat - ${describer}.gtf > ${describer}_final.gtf
done
#making log2 comparisons
bamCompare -b1 A1_lsd2ftp.sorted.bam -b2 C1_wtbackground.sorted.bam --operation log2 -p 10 -bs 1 --scaleFactorsMethod readCount --pseudocount 1 -of bigwig -o A1_lsd2ftp_C1_wtbackground_log2_combined.bw
bamCompare -b1 A2_lsd2deltacftp.sorted.bam -b2 C1_wtbackground.sorted.bam --operation log2 -p 10 -bs 1 --scaleFactorsMethod readCount --pseudocount 1 -of bigwig -o A2_lsd2deltacftp_C1_wtbackground_log2_combined.bw
bamCompare -b1 A3_lsd2ftpclr4del.sorted.bam -b2 C1_wtbackground.sorted.bam --operation log2 -p 10 -bs 1 --scaleFactorsMethod readCount --pseudocount 1 -of bigwig -o A3_lsd2ftpclr4del_C1_wtbackground_log2_combined.bw
bamCompare -b1 A4_lsd2ftpset1del.sorted.bam -b2 C1_wtbackground.sorted.bam --operation log2 -p 10 -bs 1 --scaleFactorsMethod readCount --pseudocount 1 -of bigwig -o A4_lsd2ftpset1del_C1_wtbackground_log2_combined.bw
bamCompare -b1 A5_lsd2ftpclr4set1del.sorted.bam -b2 C1_wtbackground.sorted.bam --operation log2 -p 10 -bs 1 --scaleFactorsMethod readCount --pseudocount 1 -of bigwig -o A5_lsd2ftpclr4set1del_C1_wtbackground_log2_combined.bw
bamCompare -b1 C2_lsd1ftp.sorted.bam -b2 C1_wtbackground.sorted.bam --operation log2 -p 10 -bs 1 --scaleFactorsMethod readCount --pseudocount 1 -of bigwig -o C2_lsd1ftp_C1_wtbackground_log2_combined.bw
bamCompare -b1 C3_lsd1deltahmgftp.sorted.bam -b2 C1_wtbackground.sorted.bam --operation log2 -p 10 -bs 1 --scaleFactorsMethod readCount --pseudocount 1 -of bigwig -o C3_lsd1deltahmgftp_C1_wtbackground_log2_combined.bw
bamCompare -b1 C4_lsd1ftpclr4del.sorted.bam -b2 C1_wtbackground.sorted.bam --operation log2 -p 10 -bs 1 --scaleFactorsMethod readCount --pseudocount 1 -of bigwig -o C4_lsd1ftpclr4del_C1_wtbackground_log2_combined.bw
bamCompare -b1 C5_lsd1ftpset1del.sorted.bam -b2 C1_wtbackground.sorted.bam --operation log2 -p 10 -bs 1 --scaleFactorsMethod readCount --pseudocount 1 -of bigwig -o C5_lsd1ftpset1del_C1_wtbackground_log2_combined.bw
bamCompare -b1 C6_lsd1ftpclr4set1del.sorted.bam -b2 C1_wtbackground.sorted.bam --operation log2 -p 10 -bs 1 --scaleFactorsMethod readCount --pseudocount 1 -of bigwig -o C6_lsd1ftpclr4set1del_C1_wtbackground_log2_combined.bw
#Using hclust to refine peaks of interest
for sample in *.bw
do
echo $sample
    describer=$(echo ${sample} |sed 's/_C1_wtbackground_log2_combined.bw//')
    echo $describer
    computeMatrix scale-regions -b 500 -a 500 --regionBodyLength 1500 -R ${describer}_final.gtf -S $sample -o ${describer}_scale500bp.txt.gz -p 10 -bs 1 --smartLabels --exonID gene -q
    #plotProfile -m ${describer}_scale500bp.txt.gz -out ${describer}_scale500bp.png --perGroup --plotTitle "" --regionsLabel ${describer} --legendLocation upper-right --samplesLabel "Compare to input" -y "log2" --yMin 0 --yMax 4 --dpi 300 --colors blue --numPlotsPerRow 1
    #doing this to focus on the TSS enrichment genes, putting into clusters, cluster 1 focuses on TSS enrichment
    #will output that enrichment file in bed, make compute matrix ONLY with cluster1 regions, then plotprofile of those
    plotHeatmap -m ${describer}_scale500bp.txt.gz -out ${describer}_scalehclust500bp.png --hclust 2 --dpi 300 --outFileSortedRegions ${describer}_hclust.bed
    #only sort by cluster 1 from the bed file
    awk '$13 == "cluster_1"' ${describer}_hclust.bed > ${describer}_cluster1.bed
    computeMatrix scale-regions -b 500 -a 500 --regionBodyLength 1500 -R ${describer}_cluster1.bed -S $sample -o ${describer}_scalecluster1_500bp.txt.gz -p 10 -bs 1 --smartLabels --exonID gene -q
    plotHeatmap -m ${describer}_scalecluster1_500bp.txt.gz -out ${describer}_scalehcluster1_500bp.png --dpi 300 --heatmapWidth 6 --regionsLabel ${describer} --plotTitle "" --samplesLabel "Compare to input"
done
#lsd1 vs lsd1 mutant
computeMatrix scale-regions -b 500 -a 500 --regionBodyLength 1500 -R C2_lsd1ftp_cluster1.bed -S C2_lsd1ftp_C1_wtbackground_log2_combined.bw C3_lsd1deltahmgftp_C1_wtbackground_log2_combined.bw -o lsd1comparesamples_scalecluster1_500bp_new.txt.gz -p 10 -bs 1 --smartLabels --exonID gene -q
plotProfile -m lsd1comparesamples_scalecluster1_500bp_new.txt.gz -out lsd1comparesamples_scalehcluster1_500bp_profile_new.png --dpi 300 --plotWidth 8 --samplesLabel "Lsd1" "Lsd1 mutant" --plotTitle "" --perGroup --regionsLabel "Compare to input" --yMax 6 --colors black red --yMax 6
#lsd2 vs lsd2 mutant
computeMatrix scale-regions -b 500 -a 500 --regionBodyLength 1500 -R A1_lsd2ftp_cluster1.bed -S A1_lsd2ftp_C1_wtbackground_log2_combined.bw A2_lsd2deltacftp_C1_wtbackground_log2_combined.bw -o lsd2comparesamples_scalecluster1_500bp_new.txt.gz -p 10 -bs 1 --smartLabels --exonID gene -q
plotProfile -m lsd2comparesamples_scalecluster1_500bp_new.txt.gz -out lsd2comparesamples_scalehcluster1500bp_profile_new.png --dpi 300 --plotWidth 8 --samplesLabel "Lsd2" "Lsd2 mutant" --plotTitle "" --perGroup --regionsLabel "Compare to input" --yMax 6 --colors black red --yMax 6
#lsd1 with other single mutants
computeMatrix scale-regions -b 500 -a 500 --regionBodyLength 1500 -R C2_lsd1ftp_cluster1.bed -S C2_lsd1ftp_C1_wtbackground_log2_combined.bw C3_lsd1deltahmgftp_C1_wtbackground_log2_combined.bw C4_lsd1ftpclr4del_C1_wtbackground_log2_combined.bw C5_lsd1ftpset1del_C1_wtbackground_log2_combined.bw -o lsd1compare_lsd1mutclr4delset1del_scalecluster1_500bp_new.txt.gz -p 10 -bs 1 --smartLabels --exonID gene -q
plotProfile -m lsd1compare_lsd1mutclr4delset1del_scalecluster1_500bp_new.txt.gz -out lsd1compare_lsd1mutclr4delset1del_scalecluster1_500bp_profile_new.png --dpi 300 --plotWidth 8 --samplesLabel "lsd1-ftp" "lsd1-dHMG-ftp" "lsd1-ftp-clr4d" "lsd1-ftp-set1d" --plotTitle "" --perGroup --regionsLabel "Compare to input" --yMax 6 --colors black red mediumblue purple --yMax 6
#lsd2 with other single mutants
computeMatrix scale-regions -b 500 -a 500 --regionBodyLength 1500 -R A1_lsd2ftp_cluster1.bed -S A1_lsd2ftp_C1_wtbackground_log2_combined.bw A2_lsd2deltacftp_C1_wtbackground_log2_combined.bw A5_lsd2ftpclr4del_C1_wtbackground_log2_combined.bw A4_lsd2ftpset1del_C1_wtbackground_log2_combined.bw -o lsd2compare_lsd2mutclr4delset1del_scalecluster1_500bp_new.txt.gz -p 10 -bs 1 --smartLabels --exonID gene -q
plotProfile -m lsd2compare_lsd2mutclr4delset1del_scalecluster1_500bp_new.txt.gz -out lsd2compare_lsd2mutclr4delset1del_scalecluster1_500bp_profile_new.png --dpi 300 --plotWidth 8 --samplesLabel "lsd2-ftp" "lsd2-dHMG-ftp" "lsd2-ftp-clr4d" "lsd2-ftp-set1d" --plotTitle "" --perGroup --regionsLabel "Compare to input" --yMax 6 --colors black red mediumblue purple --yMax 6