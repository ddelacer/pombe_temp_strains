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
#want the mean of each strain type
bigwigCompare -b1 A1_WT_1_cpm.bw -b2 B1_WT_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A1B1_WT_mean.bw
bigwigCompare -b1 A2_WT_37C_1_cpm.bw -b2 B2_WT_37C_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A2B2_WT_37C_mean.bw
bigwigCompare -b1 A3_lsd1_dHMG_1_cpm.bw -b2 B3_lsd1_dHMG_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A3B3_lsd1_dHMG_mean.bw
bigwigCompare -b1 A4_lsd1_dHMG_37C_1_cpm.bw -b2 B4_lsd1_dHMG_37C_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A4B4_dHMG_37C_mean.bw
bigwigCompare -b1 A5_lsd2_dC_1_cpm.bw -b2 B5_lsd2_dC_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A5B5_lsd2_dC_mean.bw
bigwigCompare -b1 A6_lsd2_dC_37C_1_cpm.bw -b2 B6_lsd2_dC_37C_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A6B6_lsd2_dC_37C_mean.bw
bigwigCompare -b1 A7_clr4d_1_cpm.bw -b2 B7_clr4d_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A7B7_clr4d_mean.bw
bigwigCompare -b1 A8_clr4d_lsd1_dHMG_1_cpm.bw -b2 B8_clr4d_lsd1_dHMG_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A8B8_clr4d_lsd1_dHMG_mean.bw
bigwigCompare -b1 A9_set1d_1_cpm.bw -b2 B9_set1d_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A9B9_set1d_mean.bw
bigwigCompare -b1 A10_set1d_37C_1_cpm.bw -b2 B10_set1d_37C_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A10B10_set1d_37C_mean.bw
bigwigCompare -b1 A11_set1d_lsd2_dC_1_cpm.bw -b2 B11_set1d_lsd2_dC_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A11B11_set1d_lsd2_dC_mean.bw
bigwigCompare -b1 A12_clr6_1_1_cpm.bw -b2 B12_clr6_1_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A12B12_clr6_mean.bw
bigwigCompare -b1 A13_clr6_1_lsd1_dHMG_1_cpm.bw -b2 B13_clr6_1_lsd1_dHMG_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A13B13_clr6_1_lsd1_mean.bw
bigwigCompare -b1 A14_clr6_1_lsd2_dC_1_cpm.bw -b2 B14_clr6_1_lsd2_dC_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A14B14_clr6_1_lsd2_mean.bw
bigwigCompare -b1 A15_sir2d_1_cpm.bw -b2 B15_sir2d_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A15B15_sir2d_mean.bw
bigwigCompare -b1 A16_sir2d_lsd1_dHMG_1_cpm.bw -b2 B16_sir2d_lsd1_dHMG_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A16B16_sir2d_lsd1_dHMG_mean.bw
bigwigCompare -b1 A17_sir2d_lsd2_dC_1_cpm.bw -b2 B17_sir2d_lsd2_dC_2_cpm.bw --operation mean -p 10 -of bigwig -bs 1 -o A17B17_sir2d_lsd2_dC_mean.bw
#forward strand
cd /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/bam_files/coverage/forward/
bigwigCompare -b1 A1_WT_1_cpm_fwd.bw -b2 B1_WT_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A1B1_WT_mean.bw
bigwigCompare -b1 A2_WT_37C_1_cpm_fwd.bw -b2 B2_WT_37C_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A2B2_WT_37C_mean.bw
bigwigCompare -b1 A3_lsd1_dHMG_1_cpm_fwd.bw -b2 B3_lsd1_dHMG_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A3B3_lsd1_dHMG_mean.bw
bigwigCompare -b1 A4_lsd1_dHMG_37C_1_cpm_fwd.bw -b2 B4_lsd1_dHMG_37C_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A4B4_dHMG_37C_mean.bw
bigwigCompare -b1 A5_lsd2_dC_1_cpm_fwd.bw -b2 B5_lsd2_dC_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A5B5_lsd2_dC_mean.bw
bigwigCompare -b1 A6_lsd2_dC_37C_1_cpm_fwd.bw -b2 B6_lsd2_dC_37C_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A6B6_lsd2_dC_37C_mean.bw
bigwigCompare -b1 A7_clr4d_1_cpm_fwd.bw -b2 B7_clr4d_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A7B7_clr4d_mean.bw
bigwigCompare -b1 A8_clr4d_lsd1_dHMG_1_cpm_fwd.bw -b2 B8_clr4d_lsd1_dHMG_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A8B8_clr4d_lsd1_dHMG_mean.bw
bigwigCompare -b1 A9_set1d_1_cpm_fwd.bw -b2 B9_set1d_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A9B9_set1d_mean.bw
bigwigCompare -b1 A10_set1d_37C_1_cpm_fwd.bw -b2 B10_set1d_37C_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A10B10_set1d_37C_mean.bw
bigwigCompare -b1 A11_set1d_lsd2_dC_1_cpm_fwd.bw -b2 B11_set1d_lsd2_dC_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A11B11_set1d_lsd2_dC_mean.bw
bigwigCompare -b1 A12_clr6_1_1_cpm_fwd.bw -b2 B12_clr6_1_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A12B12_clr6_mean.bw
bigwigCompare -b1 A13_clr6_1_lsd1_dHMG_1_cpm_fwd.bw -b2 B13_clr6_1_lsd1_dHMG_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A13B13_clr6_1_lsd1_mean.bw
bigwigCompare -b1 A14_clr6_1_lsd2_dC_1_cpm_fwd.bw -b2 B14_clr6_1_lsd2_dC_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A14B14_clr6_1_lsd2_mean.bw
bigwigCompare -b1 A15_sir2d_1_cpm_fwd.bw -b2 B15_sir2d_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A15B15_sir2d_mean.bw
bigwigCompare -b1 A16_sir2d_lsd1_dHMG_1_cpm_fwd.bw -b2 B16_sir2d_lsd1_dHMG_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A16B16_sir2d_lsd1_dHMG_mean.bw
bigwigCompare -b1 A17_sir2d_lsd2_dC_1_cpm_fwd.bw -b2 B17_sir2d_lsd2_dC_2_cpm_fwd.bw --operation mean -p 10 -of bigwig -bs 1 -o A17B17_sir2d_lsd2_dC_mean.bw
mv *mean.bw mean_coverage
cd /deac/bio/peaseGrp/ddelacer/pombe/rna_seq_7.27.2021/usftp21.novogene.com/bam_files/coverage/reverse/
bigwigCompare -b1 A1_WT_1_cpm_rev.bw -b2 B1_WT_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A1B1_WT_mean.bw
bigwigCompare -b1 A2_WT_37C_1_cpm_rev.bw -b2 B2_WT_37C_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A2B2_WT_37C_mean.bw
bigwigCompare -b1 A3_lsd1_dHMG_1_cpm_rev.bw -b2 B3_lsd1_dHMG_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A3B3_lsd1_dHMG_mean.bw
bigwigCompare -b1 A4_lsd1_dHMG_37C_1_cpm_rev.bw -b2 B4_lsd1_dHMG_37C_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A4B4_dHMG_37C_mean.bw
bigwigCompare -b1 A5_lsd2_dC_1_cpm_rev.bw -b2 B5_lsd2_dC_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A5B5_lsd2_dC_mean.bw
bigwigCompare -b1 A6_lsd2_dC_37C_1_cpm_rev.bw -b2 B6_lsd2_dC_37C_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A6B6_lsd2_dC_37C_mean.bw
bigwigCompare -b1 A7_clr4d_1_cpm_rev.bw -b2 B7_clr4d_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A7B7_clr4d_mean.bw
bigwigCompare -b1 A8_clr4d_lsd1_dHMG_1_cpm_rev.bw -b2 B8_clr4d_lsd1_dHMG_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A8B8_clr4d_lsd1_dHMG_mean.bw
bigwigCompare -b1 A9_set1d_1_cpm_rev.bw -b2 B9_set1d_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A9B9_set1d_mean.bw
bigwigCompare -b1 A10_set1d_37C_1_cpm_rev.bw -b2 B10_set1d_37C_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A10B10_set1d_37C_mean.bw
bigwigCompare -b1 A11_set1d_lsd2_dC_1_cpm_rev.bw -b2 B11_set1d_lsd2_dC_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A11B11_set1d_lsd2_dC_mean.bw
bigwigCompare -b1 A12_clr6_1_1_cpm_rev.bw -b2 B12_clr6_1_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A12B12_clr6_mean.bw
bigwigCompare -b1 A13_clr6_1_lsd1_dHMG_1_cpm_rev.bw -b2 B13_clr6_1_lsd1_dHMG_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A13B13_clr6_1_lsd1_mean.bw
bigwigCompare -b1 A14_clr6_1_lsd2_dC_1_cpm_rev.bw -b2 B14_clr6_1_lsd2_dC_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A14B14_clr6_1_lsd2_mean.bw
bigwigCompare -b1 A15_sir2d_1_cpm_rev.bw -b2 B15_sir2d_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A15B15_sir2d_mean.bw
bigwigCompare -b1 A16_sir2d_lsd1_dHMG_1_cpm_rev.bw -b2 B16_sir2d_lsd1_dHMG_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A16B16_sir2d_lsd1_dHMG_mean.bw
bigwigCompare -b1 A17_sir2d_lsd2_dC_1_cpm_rev.bw -b2 B17_sir2d_lsd2_dC_2_cpm_rev.bw --operation mean -p 10 -of bigwig -bs 1 -o A17B17_sir2d_lsd2_dC_mean.bw
mv *mean.bw mean_coverage