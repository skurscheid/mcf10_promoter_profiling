#!/bin/bash

export cli_params="--missingDataAsZero --skipZeros --smartLabels --referencePoint center --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --sortRegions keep --numberOfProcessors 32"
export scoreFiles="deeptools/bigwigCompare/MCF10A/H3K4me1_vs_Input.bw\
                deeptools/bigwigCompare/MCF10A/H3K4me3_vs_Input.bw\
                deeptools/bigwigCompare/MCF10A/H3K79me2_vs_Input.bw\
                deeptools/bigwigCompare/MCF10A/H3K27me3_vs_Input.bw\
                deeptools/bigwigCompare/MCF10A/H3K27ac_vs_Input.bw\
                deeptools/bigwigCompare/MCF10A/H2BK120ub1_vs_Input.bw\
                deeptools/bigwigCompare/MCF10A/H3K9me3_vs_Input.bw\
                deeptools/bigwigCompare/MCF10A/H3K9ac_vs_Input.bw\
                deeptools/bigwigCompare/MCF10A/H3K36me3_vs_Input.bw\
                deeptools/bigwigCompare/MCF10A/H3K23ac_vs_Input.bw\
                deeptools/bigwigCompare/MCF10A/H4K8ac_vs_Input.bw"

computeMatrix reference-point --regionsFileName /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1A_cluster_7.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1A_cluster_6.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1A_cluster_5.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1A_cluster_4.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1A_cluster_3.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1A_cluster_2.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1A_cluster_1.bed\
                              --scoreFileName $scoreFiles\
                              --outFileName deeptools/computeMatrix_referencepoint/MCF10A/Fig1A_matrix.gz
                              2>logs/deeptools_computeMatrix/MCF10A/Fig1A_matrix.log

plotHeatmap -m deeptools/computeMatrix_referencepoint/MCF10A/Fig1A_matrix.gz\
            -o deeptools/plotHeatmap/MCF10A/Fig1A/Fig1A_heatmap.pdf\
            --regionsLabel 7 6 5 4 3 2 1\
            --sortRegions keep\
            --whatToShow 'heatmap and colorbar'\
            --plotTitle 'Figure 1A - histone marks'\
            2>logs/deeptools_plotHeatmap/MCF10A/Fig1A_heatmap.log

plotHeatmap -m deeptools/computeMatrix_referencepoint/MCF10A/Fig1A_matrix.gz\
            -o deeptools/plotHeatmap/MCF10A/Fig1A/Fig1A_heatmap.png\
            --regionsLabel 7 6 5 4 3 2 1\
            --sortRegions keep\
            --whatToShow 'heatmap and colorbar'\
            --plotTitle 'Figure 1A - histone marks'\
            2>>logs/deeptools_plotHeatmap/MCF10A/Fig1A_heatmap.log


# Figure 1B histone marks
computeMatrix reference-point --regionsFileName /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1B_cluster_7.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1B_cluster_6.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1B_cluster_5.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1B_cluster_4.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1B_cluster_3.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1B_cluster_2.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1B_cluster_1.bed\
                              --scoreFileName $scoreFiles\
                              $cli_params\
                              --outFileName deeptools/computeMatrix_referencepoint/MCF10A/Fig1B_matrix.gz\
                              2>logs/deeptools_computeMatrix/MCF10A/Fig1B_matrix.log

plotHeatmap -m deeptools/computeMatrix_referencepoint/MCF10A/Fig1B_matrix.gz\
            -o deeptools/plotHeatmap/MCF10A/Fig1B/Fig1B_heatmap.pdf\
            --regionsLabel 7 6 5 4 3 2 1\
            --sortRegions keep\
            --plotTitle 'Figure 1B - histone marks'\
            2>logs/deeptools_plotHeatmap/MCF10A/Fig1B_heatmap.log

plotHeatmap -m deeptools/computeMatrix_referencepoint/MCF10A/Fig1B_matrix.gz\
            -o deeptools/plotHeatmap/MCF10A/Fig1B/Fig1B_heatmap.png\
            --regionsLabel 7 6 5 4 3 2 1\
            --sortRegions keep\
            --plotTitle 'Figure 1B - histone marks'\
            2>>logs/deeptools_plotHeatmap/MCF10A/Fig1B_heatmap.log

# Figure 1C histone marks
computeMatrix reference-point --regionsFileName /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1C_cluster_7.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1C_cluster_6.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1C_cluster_5.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1C_cluster_4.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1C_cluster_3.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1C_cluster_2.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1C_cluster_1.bed\
                              --scoreFileName $scoreFiles\
                              --outFileName deeptools/computeMatrix_referencepoint/MCF10A/Fig1C_matrix.gz\
                              $cli_params\
                              2>logs/deeptools_computeMatrix/MCF10A/Fig1C_matrix.log

plotHeatmap -m deeptools/computeMatrix_referencepoint/MCF10A/Fig1C_matrix.gz\
            -o deeptools/plotHeatmap/MCF10A/Fig1C/Fig1C_heatmap.pdf\
            --regionsLabel 7 6 5 4 3 2 1\
            --sortRegions keep\
            --plotTitle 'Figure 1C - histone marks'\
            2>logs/deeptools_plotHeatmap/MCF10A/Fig1C_heatmap.log

plotHeatmap -m deeptools/computeMatrix_referencepoint/MCF10A/Fig1C_matrix.gz\
            -o deeptools/plotHeatmap/MCF10A/Fig1C/Fig1C_heatmap.png\
            --regionsLabel 7 6 5 4 3 2 1\
            --sortRegions keep\
            --plotTitle 'Figure 1C - histone marks'\
            2>>logs/deeptools_plotHeatmap/MCF10A/Fig1C_heatmap.log

# Figure 1D histone marks
computeMatrix reference-point --regionsFileName /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1D_cluster_7.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1D_cluster_6.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1D_cluster_5.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1D_cluster_4.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1D_cluster_3.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1D_cluster_2.bed\
                                                /home/150/sxk150/data/PromoterSeqCap/annotations/Fig1D_cluster_1.bed\
                              --scoreFileName $scoreFiles\
                              --outFileName deeptools/computeMatrix_referencepoint/MCF10A/Fig1D_matrix.gz\
                              $cli_params\
                              2>logs/deeptools_computeMatrix/MCF10A/Fig1D_matrix.log
        
plotHeatmap -m deeptools/computeMatrix_referencepoint/MCF10A/Fig1D_matrix.gz\
            -o deeptools/plotHeatmap/MCF10A/Fig1D/Fig1D_heatmap.pdf\
            --regionsLabel 7 6 5 4 3 2 1\
            --sortRegions keep\
            --plotTitle 'Figure 1D - histone marks'\
            2>logs/deeptools_plotHeatmap/MCF10A/Fig1D_heatmap.log

plotHeatmap -m deeptools/computeMatrix_referencepoint/MCF10A/Fig1D_matrix.gz\
            -o deeptools/plotHeatmap/MCF10A/Fig1D/Fig1D_heatmap.png\
            --regionsLabel 7 6 5 4 3 2 1\
            --sortRegions keep\
            --plotTitle 'Figure 1D - histone marks'\
            2>>logs/deeptools_plotHeatmap/MCF10A/Fig1D_heatmap.log