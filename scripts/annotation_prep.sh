#!/bin/bash

# preparing single files per cluster for plotting of the profile/heatmaps with deeptools
grep -v sort1 Fig1C_sorting.tsv | cut -f 1,2,4 | sed s/\;/\\t/g | cut -f 1,2,3,5,7,8,4 | awk 'BEGIN {OFS="\t";}{$7="cluster_"$7; print $0;}' | sort -k 7 -k 6V - > Fig1C_sorting_group_sorted.bed
for i in $(cut -f 7 Fig1C_sorting_group_sorted.bed|uniq); do grep $i Fig1C_sorting_group_sorted.bed > Fig1D_${i}.bed; done

grep -v sort1 Fig1D_sorting.tsv | cut -f 1,2,4 | sed s/\;/\\t/g | cut -f 1,2,3,5,7,8,4 | awk 'BEGIN {OFS="\t";}{$7="cluster_"$7; print $0;}' | sort -k 7 -k 6V - > Fig1D_sorting_group_sorted.bed
for i in $(cut -f 7 Fig1D_sorting_group_sorted.bed|uniq); do grep $i Fig1D_sorting_group_sorted.bed > Fig1D_${i}.bed; done

grep -v sort1 Fig1A_Total_10A_Input_k7sorting.tsv | cut -f 1,2,4 | sed s/\;/\\t/g | cut -f 1,2,3,5,7,8,4 | awk 'BEGIN {OFS="\t";}{$7="cluster_"$7; print $0;}' | sort -k 7 -k 6V - > Fig1A_Total_10A_Input_k7sorting.bed
for i in $(cut -f 7 Fig1A_Total_10A_Input_k7sorting.bed|uniq); do grep $i Fig1A_Total_10A_Input_k7sorting.bed > Fig1A_${i}.bed; done

grep -v sort1 Fig1B_TOTAL_10A_H2AZ_k7sorting.tsv | cut -f 1,2,4 | sed s/\;/\\t/g | cut -f 1,2,3,5,7,8,4 | awk 'BEGIN {OFS="\t";}{$7="cluster_"$7; print $0;}' | sort -k 7 -k 6V - > Fig1B_TOTAL_10A_H2AZ_k7sorting.bed
for i in $(cut -f 7 Fig1B_TOTAL_10A_H2AZ_k7sorting.bed|uniq); do grep $i Fig1B_TOTAL_10A_H2AZ_k7sorting.bed > Fig1B_${i}.bed; done